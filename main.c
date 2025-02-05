#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#define BUFFER_SIZE (1024 * 1024)  // 1MB buffer
#define OVERLAP_SIZE 1000          // Size of overlap between chunks
#define MAX_PATTERNS 1000
#define MAX_PATTERN_LENGTH 1000
#define MAX_HEADER_LENGTH 1000
#define ALPHABET_SIZE 256

typedef struct {
    char name[256];
    char sequence[MAX_PATTERN_LENGTH];
    char preprocessed[MAX_PATTERN_LENGTH];
    size_t length;
} Pattern;

typedef struct {
    char header[MAX_HEADER_LENGTH];
    char* sequence;
    size_t position;
    size_t global_position;  // Position in original sequence
    char pattern_name[256];
    char pattern_sequence[MAX_PATTERN_LENGTH];
    int strand;  // 0 for forward, 1 for reverse
} FastaMatch;

typedef struct {
    char* data;
    size_t length;
    size_t capacity;
    size_t global_offset;  // Track position in original sequence
} SequenceBuffer;

typedef struct {
    SequenceBuffer forward;
    SequenceBuffer reverse;
    char* preprocessed;
    char* rev_preprocessed;
    size_t max_pattern_length;
} ChunkProcessor;

typedef struct ACNode {
    struct ACNode* children[ALPHABET_SIZE];
    struct ACNode* failure;
    int* pattern_indices;
    int pattern_count;
    int capacity;
} ACNode;

typedef struct {
    ACNode** items;
    int front;
    int rear;
    int capacity;
} Queue;

Queue* create_queue(int size) {
    Queue* queue = malloc(sizeof(Queue));
    queue->items = malloc(size * sizeof(ACNode*));
    queue->front = queue->rear = -1;
    queue->capacity = size;
    return queue;
}

void enqueue(Queue* queue, ACNode* item) {
    if (queue->rear == -1) {
        queue->front = queue->rear = 0;
    } else {
        queue->rear = (queue->rear + 1) % queue->capacity;
    }
    queue->items[queue->rear] = item;
}

ACNode* dequeue(Queue* queue) {
    ACNode* item = queue->items[queue->front];
    if (queue->front == queue->rear) {
        queue->front = queue->rear = -1;
    } else {
        queue->front = (queue->front + 1) % queue->capacity;
    }
    return item;
}

int is_empty(Queue* queue) {
    return queue->front == -1;
}

void free_queue(Queue* queue) {
    free(queue->items);
    free(queue);
}

ACNode* create_node() {
    ACNode* node = calloc(1, sizeof(ACNode));
    node->pattern_indices = malloc(10 * sizeof(int));
    node->capacity = 10;
    node->pattern_count = 0;
    return node;
}

void add_pattern_to_node(ACNode* node, int pattern_index) {
    if (node->pattern_count == node->capacity) {
        node->capacity *= 2;
        node->pattern_indices = realloc(node->pattern_indices, node->capacity * sizeof(int));
    }
    node->pattern_indices[node->pattern_count++] = pattern_index;
}

void free_node(ACNode* node) {
    if (!node) return;
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        free_node(node->children[i]);
    }
    free(node->pattern_indices);
    free(node);
}

ACNode* build_automaton(Pattern* patterns, int pattern_count, int ignore_case) {
    ACNode* root = create_node();
    
    // Build trie
    for (int i = 0; i < pattern_count; i++) {
        ACNode* current = root;
        char* pattern = patterns[i].preprocessed;
        
        for (size_t j = 0; pattern[j]; j++) {
            unsigned char c = pattern[j];
            if (!current->children[c]) {
                current->children[c] = create_node();
            }
            current = current->children[c];
        }
        add_pattern_to_node(current, i);
    }
    
    Queue* queue = create_queue(pattern_count * MAX_PATTERN_LENGTH);
    
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (root->children[i]) {
            root->children[i]->failure = root;
            enqueue(queue, root->children[i]);
        }
    }
    
    while (!is_empty(queue)) {
        ACNode* current = dequeue(queue);
        
        for (int i = 0; i < ALPHABET_SIZE; i++) {
            if (current->children[i]) {
                ACNode* state = current->failure;
                while (state && !state->children[i]) {
                    state = state->failure;
                }
                
                current->children[i]->failure = state ? state->children[i] : root;
                
                for (int j = 0; j < current->children[i]->failure->pattern_count; j++) {
                    add_pattern_to_node(current->children[i], 
                                      current->children[i]->failure->pattern_indices[j]);
                }
                
                enqueue(queue, current->children[i]);
            }
        }
    }
    
    free_queue(queue);
    return root;
}

int read_patterns(const char* filepath, Pattern* patterns, int ignore_case) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        perror("Error opening patterns file");
        return -1;
    }
    
    char line[MAX_PATTERN_LENGTH];
    int count = 0;
    
    // Skip header
    fgets(line, sizeof(line), file);
    
    while (fgets(line, sizeof(line), file) && count < MAX_PATTERNS) {
        char* name = strtok(line, ",");
        char* seq = strtok(NULL, "\n");
        
        if (name && seq) {
            strncpy(patterns[count].name, name, sizeof(patterns[count].name) - 1);
            strncpy(patterns[count].sequence, seq, sizeof(patterns[count].sequence) - 1);
            patterns[count].length = strlen(seq);
            
            for (size_t i = 0; i < patterns[count].length; i++) {
                patterns[count].preprocessed[i] = ignore_case ? 
                    tolower(seq[i]) : seq[i];
            }
            patterns[count].preprocessed[patterns[count].length] = '\0';
            
            count++;
        }
    }
    
    fclose(file);
    return count;
}

char complement(char c) {
    switch(toupper(c)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default: return 'N';
    }
}

void init_sequence_buffer(SequenceBuffer* buf, size_t initial_capacity) {
    buf->data = malloc(initial_capacity);
    buf->length = 0;
    buf->capacity = initial_capacity;
    buf->global_offset = 0;
}

ChunkProcessor* create_chunk_processor(size_t chunk_size, size_t max_pattern_len) {
    ChunkProcessor* processor = malloc(sizeof(ChunkProcessor));
    init_sequence_buffer(&processor->forward, chunk_size + max_pattern_len);
    init_sequence_buffer(&processor->reverse, chunk_size + max_pattern_len);
    processor->preprocessed = malloc(chunk_size + max_pattern_len);
    processor->rev_preprocessed = malloc(chunk_size + max_pattern_len);
    processor->max_pattern_length = max_pattern_len;
    return processor;
}

void process_chunk(ChunkProcessor* processor, const char* new_data, size_t new_len,
                  ACNode* root, Pattern* patterns, const char* header,
                  size_t context, int ignore_case, FastaMatch** matches, size_t* match_count,
                  int is_final_chunk) {
    
    if (processor->forward.length + new_len > processor->forward.capacity) {
        size_t new_cap = processor->forward.capacity * 2;
        processor->forward.data = realloc(processor->forward.data, new_cap);
        processor->preprocessed = realloc(processor->preprocessed, new_cap);
        processor->forward.capacity = new_cap;
    }
    
    memcpy(processor->forward.data + processor->forward.length, new_data, new_len);
    
    for (size_t i = 0; i < new_len; i++) {
        processor->preprocessed[processor->forward.length + i] = 
            ignore_case ? tolower(new_data[i]) : new_data[i];
    }
    
    processor->forward.length += new_len;
    
    if (processor->forward.length >= OVERLAP_SIZE || is_final_chunk) {
        size_t process_len = is_final_chunk ? processor->forward.length : 
                            processor->forward.length - processor->max_pattern_length;
        
        ACNode* current = root;
        for (size_t i = 0; i < process_len; i++) {
            unsigned char c = processor->preprocessed[i];
            while (current != root && !current->children[c]) {
                current = current->failure;
            }
            current = current->children[c] ? current->children[c] : root;
            
            if (current->pattern_count > 0) {
                for (int j = 0; j < current->pattern_count; j++) {
                    int pattern_idx = current->pattern_indices[j];
                    size_t global_pos = processor->forward.global_offset + i;
                    
                    *matches = realloc(*matches, (*match_count + 1) * sizeof(FastaMatch));
                    FastaMatch* match = &(*matches)[*match_count];
                    
                    strncpy(match->header, header, MAX_HEADER_LENGTH - 1);
                    match->position = i - patterns[pattern_idx].length + 1;
                    match->global_position = global_pos - patterns[pattern_idx].length + 1;
                    strncpy(match->pattern_name, patterns[pattern_idx].name, 255);
                    strncpy(match->pattern_sequence, patterns[pattern_idx].sequence, MAX_PATTERN_LENGTH - 1);
                    
                    size_t ctx_start = match->position > context ? match->position - context : 0;
                    size_t ctx_end = (match->position + patterns[pattern_idx].length + context) < process_len ? 
                                    match->position + patterns[pattern_idx].length + context : process_len;
                    
                    size_t ctx_len = ctx_end - ctx_start;
                    match->sequence = malloc(ctx_len + 1);
                    strncpy(match->sequence, &processor->forward.data[ctx_start], ctx_len);
                    match->sequence[ctx_len] = '\0';
                    match->strand = 0;
                    
                    (*match_count)++;
                }
            }
        }
        
        // Process reverse complement
        char* rev_seq = malloc(process_len);
        for (size_t i = 0; i < process_len; i++) {
            rev_seq[i] = complement(processor->forward.data[process_len - 1 - i]);
            processor->rev_preprocessed[i] = ignore_case ? tolower(rev_seq[i]) : rev_seq[i];
        }
        
        current = root;
        for (size_t i = 0; i < process_len; i++) {
            unsigned char c = processor->rev_preprocessed[i];
            while (current != root && !current->children[c]) {
                current = current->failure;
            }
            current = current->children[c] ? current->children[c] : root;
            
            if (current->pattern_count > 0) {
                for (int j = 0; j < current->pattern_count; j++) {
                    int pattern_idx = current->pattern_indices[j];
                    size_t global_pos = processor->forward.global_offset + process_len - 1 - i;
                    
                    *matches = realloc(*matches, (*match_count + 1) * sizeof(FastaMatch));
                    FastaMatch* match = &(*matches)[*match_count];
                    
                    strncpy(match->header, header, MAX_HEADER_LENGTH - 1);
                    match->position = global_pos;
                    match->global_position = global_pos;
                    strncpy(match->pattern_name, patterns[pattern_idx].name, 255);
                    strncpy(match->pattern_sequence, patterns[pattern_idx].sequence, MAX_PATTERN_LENGTH - 1);
                    
                    size_t rev_ctx_start = i > context ? i - context : 0;
                    size_t rev_ctx_end = i + patterns[pattern_idx].length + context;
                    if (rev_ctx_end > process_len) rev_ctx_end = process_len;
                    
                    size_t ctx_len = rev_ctx_end - rev_ctx_start;
                    match->sequence = malloc(ctx_len + 1);
                    for (size_t k = 0; k < ctx_len; k++) {
                        match->sequence[k] = rev_seq[rev_ctx_start + k];
                    }
                    match->sequence[ctx_len] = '\0';
                    match->strand = 1;
                    
                    (*match_count)++;
                }
            }
        }
        
        free(rev_seq);
        
        if (!is_final_chunk) {
            memmove(processor->forward.data, 
                   processor->forward.data + process_len,
                   processor->forward.length - process_len);
            processor->forward.length -= process_len;
            processor->forward.global_offset += process_len;
        }
    }
}

void free_chunk_processor(ChunkProcessor* processor) {
    free(processor->forward.data);
    free(processor->reverse.data);
    free(processor->preprocessed);
    free(processor->rev_preprocessed);
    free(processor);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <fasta_file> <patterns_file> [context] [sequence_only] [ignore_case]\n", 
                argv[0]);
        return 1;
    }
    
    const char* fasta_file = argv[1];
    const char* patterns_file = argv[2];
    size_t context = argc > 3 ? atoi(argv[3]) : 0;
    int sequence_only = argc > 4 ? atoi(argv[4]) : 0;
    int ignore_case = argc > 5 ? atoi(argv[5]) : 0;
    
    Pattern patterns[MAX_PATTERNS];
    int pattern_count = read_patterns(patterns_file, patterns, ignore_case);
    if (pattern_count < 0) return 1;
    
    fprintf(stderr,"Loaded %d patterns\n", pattern_count);
    
    size_t max_pattern_len = 0;
    for (int i = 0; i < pattern_count; i++) {
        if (patterns[i].length > max_pattern_len) {
            max_pattern_len = patterns[i].length;
        }
    }
    
    ACNode* root = build_automaton(patterns, pattern_count, ignore_case);
    
    gzFile fp = gzopen(fasta_file, "r");
    if (!fp) {
        perror("Error opening FASTA file");
        return 1;
    }
    
    ChunkProcessor* processor = create_chunk_processor(BUFFER_SIZE, max_pattern_len);
    FastaMatch* matches = NULL;
    size_t match_count = 0;
    char current_header[MAX_HEADER_LENGTH] = "";
    char buffer[BUFFER_SIZE];
    int in_header = 0;
    
    // Print CSV header
    printf("header,pattern_name,pattern_sequence,position,strand,context\n");
    
    while (1) {
        size_t bytes_read = gzread(fp, buffer, BUFFER_SIZE);
        if (bytes_read == 0) break;
        
        char* line_start = buffer;
        char* current = buffer;
        char* end = buffer + bytes_read;
        
        while (current < end) {
            if (*current == '>' || current == end - 1) {
                if (in_header) {
                    // Store header
                    size_t header_len = current - line_start;
                    if (header_len > MAX_HEADER_LENGTH - 1) 
                        header_len = MAX_HEADER_LENGTH - 1;
                    strncpy(current_header, line_start + 1, header_len);
                    current_header[header_len] = '\0';
                    current_header[strcspn(current_header, "\n")] = 0;
                } else if (current - line_start > 0) {
                    // Process sequence chunk
                    process_chunk(processor, line_start, current - line_start,
                                root, patterns, current_header, context,
                                ignore_case, &matches, &match_count,
                                bytes_read < BUFFER_SIZE);
                }
                
                in_header = (*current == '>');
                line_start = current;
            }
            current++;
        }
    }
    
    // Print matches in CSV format
    for (size_t i = 0; i < match_count; i++) {
        // Escape any commas in the fields
        char* header = strdup(matches[i].header);
        char* pattern_name = strdup(matches[i].pattern_name);
        char* pattern_sequence = strdup(matches[i].pattern_sequence);
        char* sequence = strdup(matches[i].sequence);
        
        for (char* p = header; *p; p++) if (*p == ',') *p = ';';
        for (char* p = pattern_name; *p; p++) if (*p == ',') *p = ';';
        for (char* p = pattern_sequence; *p; p++) if (*p == ',') *p = ';';
        for (char* p = sequence; *p; p++) if (*p == ',') *p = ';';
        
        printf("%s,%s,%s,%zu,%s,%s\n",
               header,
               pattern_name,
               pattern_sequence,
               matches[i].global_position,
               matches[i].strand ? "reverse" : "forward",
               sequence);
               
        free(header);
        free(pattern_name);
        free(pattern_sequence);
        free(sequence);
        free(matches[i].sequence);
    }
    
    // Cleanup
    free_chunk_processor(processor);
    free(matches);
    free_node(root);
    gzclose(fp);
    
    return 0;
}
