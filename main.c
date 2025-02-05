#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#define BUFFER_SIZE (1024 * 1024)  // 1MB buffer
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
    char pattern_name[256];
    char pattern_sequence[MAX_PATTERN_LENGTH];
    int strand;  // 0 for forward, 1 for reverse
} FastaMatch;

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

// Queue operations
Queue* create_queue(int size) {
    Queue* queue = (Queue*)malloc(sizeof(Queue));
    queue->items = (ACNode**)malloc(size * sizeof(ACNode*));
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

// Aho-Corasick node operations
ACNode* create_node() {
    ACNode* node = (ACNode*)calloc(1, sizeof(ACNode));
    node->pattern_indices = (int*)malloc(10 * sizeof(int));
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

// Build Aho-Corasick automaton
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
    
    // Build failure links using BFS
    Queue* queue = create_queue(pattern_count * MAX_PATTERN_LENGTH);
    
    // Set root's children failure to root
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
                
                // Copy pattern indices from failure state
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

// Read patterns from CSV
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
            
            // Create preprocessed version
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

// Search sequence using Aho-Corasick
char complement(char c) {
    switch(toupper(c)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';
        default: return 'N';
    }
}

void create_reverse_complement(const char* sequence, char* reverse, size_t length) {
    for (size_t i = 0; i < length; i++) {
        reverse[i] = complement(sequence[length - 1 - i]);
    }
    reverse[length] = '\0';
}

void search_sequence(const char* sequence, const char* original_sequence, size_t seq_length,
                    ACNode* root, Pattern* patterns, const char* header,
                    size_t context, int ignore_case, FastaMatch** matches, size_t* match_count) {
    ACNode* current = root;
    
    for (size_t i = 0; i < seq_length; i++) {
        unsigned char c = ignore_case ? tolower(sequence[i]) : sequence[i];
        
        while (current != root && !current->children[c]) {
            current = current->failure;
        }
        
        current = current->children[c] ? current->children[c] : root;
        
        if (current->pattern_count > 0) {
            for (int j = 0; j < current->pattern_count; j++) {
                int pattern_idx = current->pattern_indices[j];
                size_t pattern_len = patterns[pattern_idx].length;
                size_t start_pos = i - pattern_len + 1;
                
                // Extend matches array
                *matches = realloc(*matches, (*match_count + 1) * sizeof(FastaMatch));
                FastaMatch* match = &(*matches)[*match_count];
                
                strncpy(match->header, header, MAX_HEADER_LENGTH - 1);
                match->position = start_pos;
                strncpy(match->pattern_name, patterns[pattern_idx].name, 255);
                strncpy(match->pattern_sequence, patterns[pattern_idx].sequence, MAX_PATTERN_LENGTH - 1);
                
                // Extract context
                size_t ctx_start = start_pos > context ? start_pos - context : 0;
                size_t ctx_end = (start_pos + pattern_len + context) < seq_length ? 
                                 start_pos + pattern_len + context : seq_length;
                
                size_t ctx_len = ctx_end - ctx_start;
                match->sequence = malloc(ctx_len + 1);
                strncpy(match->sequence, &original_sequence[ctx_start], ctx_len);
                match->sequence[ctx_len] = '\0';
                
                match->strand = 0;  // forward strand
                (*match_count)++;
            }
        }
    }
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
    
    // Read patterns
    Pattern patterns[MAX_PATTERNS];
    int pattern_count = read_patterns(patterns_file, patterns, ignore_case);
    if (pattern_count < 0) return 1;
    
    printf("Loaded %d patterns\n", pattern_count);
    
    // Build Aho-Corasick automaton
    ACNode* root = build_automaton(patterns, pattern_count, ignore_case);
    
    // Process FASTA file
    gzFile fp = gzopen(fasta_file, "r");
    if (!fp) {
        perror("Error opening FASTA file");
        return 1;
    }
    
    char* buffer = malloc(BUFFER_SIZE);
    char* sequence = NULL;
    char* preprocessed = NULL;
    size_t seq_capacity = 0;
    size_t seq_length = 0;
    char current_header[MAX_HEADER_LENGTH] = "";
    FastaMatch* matches = NULL;
    size_t match_count = 0;
    size_t last_match_count = 0;
    
    while (gzgets(fp, buffer, BUFFER_SIZE)) {
        if (buffer[0] == '>') {
            // Process previous sequence if exists
            if (seq_length > 0) {
                sequence[seq_length] = '\0';
                preprocessed[seq_length] = '\0';
                
                // Create reverse complement
                char* reverse_seq = malloc(seq_length + 1);
                char* reverse_preprocessed = malloc(seq_length + 1);
                create_reverse_complement(sequence, reverse_seq, seq_length);
                for (size_t i = 0; i < seq_length; i++) {
                    reverse_preprocessed[i] = ignore_case ? tolower(reverse_seq[i]) : reverse_seq[i];
                }
                reverse_preprocessed[seq_length] = '\0';
                
                // Search forward strand
                last_match_count = match_count;
                search_sequence(preprocessed, sequence, seq_length,
                              root, patterns, current_header,
                              context, ignore_case, &matches, &match_count);
                
                seq_length = 0;
            }
            
            // Store new header
            strncpy(current_header, buffer + 1, MAX_HEADER_LENGTH - 1);
            current_header[strcspn(current_header, "\n")] = 0;
        } else {
            size_t line_len = strlen(buffer);
            buffer[strcspn(buffer, "\n")] = 0;
            
            // Resize buffers if needed
            if (seq_length + line_len >= seq_capacity) {
                seq_capacity = seq_capacity == 0 ? BUFFER_SIZE : seq_capacity * 2;
                sequence = realloc(sequence, seq_capacity);
                preprocessed = realloc(preprocessed, seq_capacity);
                if (!sequence || !preprocessed) {
                    fprintf(stderr, "Memory allocation failed\n");
                    return 1;
                }
            }
            
            // Add sequence data
            for (size_t i = 0; buffer[i]; i++) {
                char c = buffer[i];
                if (!isspace(c)) {
                    sequence[seq_length] = c;
                    preprocessed[seq_length] = ignore_case ? tolower(c) : c;
                    seq_length++;
                }
            }
        }
    }
    
    // Process final sequence
    if (seq_length > 0) {
        sequence[seq_length] = '\0';
        preprocessed[seq_length] = '\0';
        
        // Create reverse complement
        char* reverse_seq = malloc(seq_length + 1);
        char* reverse_preprocessed = malloc(seq_length + 1);
        create_reverse_complement(sequence, reverse_seq, seq_length);
        for (size_t i = 0; i < seq_length; i++) {
            reverse_preprocessed[i] = ignore_case ? tolower(reverse_seq[i]) : reverse_seq[i];
        }
        reverse_preprocessed[seq_length] = '\0';
        
        // Search forward strand
        search_sequence(preprocessed, sequence, seq_length,
                      root, patterns, current_header,
                      context, ignore_case, &matches, &match_count);
                      
        // Search reverse strand
        search_sequence(reverse_preprocessed, reverse_seq, seq_length,
                      root, patterns, current_header,
                      context, ignore_case, &matches, &match_count);
        
        // Mark reverse strand matches
        for (size_t i = last_match_count; i < match_count; i++) {
            matches[i].strand = 1;  // reverse strand
            // Convert position to forward strand coordinates
            matches[i].position = seq_length - matches[i].position - 1;
        }
                      
        free(reverse_seq);
        free(reverse_preprocessed);
    }
    
    // Print CSV header
    printf("header,pattern_name,pattern_sequence,position,strand,context\n");
    
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
               matches[i].position,
               matches[i].strand ? "reverse" : "forward",
               sequence);
               
        free(header);
        free(pattern_name);
        free(pattern_sequence);
        free(sequence);
        free(matches[i].sequence);
    }
    
    // Cleanup
    free(buffer);
    free(sequence);
    free(preprocessed);
    free(matches);
    free_node(root);
    gzclose(fp);
    
    return 0;
}
