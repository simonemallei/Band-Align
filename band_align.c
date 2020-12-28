#include "band_align.h"

#define INDEL '*'
#define LETTER 'A'

/* Given two int64_t, returns the one with minimum value */
int64_t min(int64_t a, int64_t b) {
	return a < b ? a : b;
}

/* Given two characters, returns the score denoted by the BLOSUM62 matrix */
int64_t get_score(char a, char b){
	int64_t index_a = 0;
	int64_t index_b = 0;
	int64_t i;
	for (i = 0; i < 24; i++){
		if (a == score_index[i])
			index_a = i;
		if (b == score_index[i])
			index_b = i;		
	}
	return matrix_score[index_a][index_b];
}

/*** Returns an upper bound of the maximum score that can be obtained using the edge of the band
   * 
   * band: size of the band
   * l2: length of s2 (longest string to align)
   * s1: shortest string to align
   */
int64_t maxim_score(int64_t band, int64_t l2, char* s1){
	int64_t l1 = strlen(s1);
	int64_t num_indel = l1 - 1 + l2 - 1;
	if (band + 1 <= l2 - 1)
		num_indel = min(num_indel, band + 1 + l1 - (l2 - band - 1));
	int64_t score = num_indel * get_score(INDEL, LETTER);
#ifdef DEBUG
	printf("Minimum number of indels: %d\n", num_indel);
#endif
	int64_t i;
	for (i = 0; i < l1; i++)
		score += get_score(s1[i], s1[i]);
	return score;
}

/*** Returns the index of (x, y) in the matrix that represents the band 
   * (-1 if (x, y) is not in the band) 
   * 
   * x: x-coordinate in the "not-band" matrix
   * y: y-coordinate in the "not-band" matrix
   * l1: length of s1, also the number of columns in the "not-band" matrix
   * band: size of the band 
   */
int64_t find_cell(int64_t x, int64_t y, int64_t l1, int64_t band){
	int64_t band_length = 2 * band + 1;
	if ((y >= x && y - x > band) || (x - y > band))
		return -1;
	if (y >= x)
		return x * band_length + y - x + band; 
	else
		return x * band_length + band - (x - y);
}

/*** Returns the best band alignment between two strings
   *
   * s1: shortest string to align
   * s2: longest string to align
   * band: size of the band 
   */
static align_s* band_align(char* s1, char* s2, int64_t band){
	int64_t l1 = strlen(s1);
    int64_t l2 = strlen(s2);
        
    cell_s* m = malloc((l2 + 1) * (band * 2 + 1) * sizeof(cell_s));
	assert(m != NULL && "Cannot allocate matrix\n");

	/* Initializing M[0, 0] in the band */
    m[find_cell(0, 0, l1, band)].cell = 0;
    m[find_cell(0, 0, l1, band)].prev_x = -1;
    m[find_cell(0, 0, l1, band)].prev_y = -1;
    m[find_cell(0, 0, l1, band)].indel_s2 = 0;
    /* Constructing M[x, 0] in the band */
	int64_t x = 1;
    for (; find_cell(x, 0, l1, band) != -1 && x <= l2; x++){
        	int64_t curr_pos = find_cell(x, 0, l1, band);
        	int64_t prev_pos = find_cell(x - 1, 0, l1, band);
    		m[curr_pos].cell = m[prev_pos].cell + get_score(INDEL, s2[x - 1]);
    		m[curr_pos].prev_x = x - 1;
        	m[curr_pos].prev_y = 0;
        	m[curr_pos].indel_s2 = m[prev_pos].indel_s2;
    }
    /* Constructing M[0, y] in the band */
	int64_t y = 1;
    for (; find_cell(0, y, l1, band) != -1 && y <= l1; y++){
    	int64_t curr_pos = find_cell(0, y, l1, band);
   		int64_t prev_pos = find_cell(0, y - 1, l1, band);
        m[curr_pos].cell = m[prev_pos].cell + get_score(s1[y - 1], INDEL);
    	m[curr_pos].prev_x = 0;
    	m[curr_pos].prev_y = y - 1;
		m[curr_pos].indel_s2 = m[prev_pos].indel_s2 + 1;
	}
	/* Constructing M[x, y] in the band */
	for (x = 1; x <= l2; x++)
        for (y = 1; y <= l1; y++){
        	int64_t curr_pos = find_cell(x, y, l1, band);
        	/* curr_pos will be -1 if M[x, y] is not in the band */
        	if (curr_pos != -1){
        		/* Uses best score from M[x-1, y-1] */
				int64_t prev_pos = find_cell(x -1, y - 1, l1, band);			
        		m[curr_pos].cell = m[prev_pos].cell + get_score(s1[y - 1], s2[x - 1]);
        		m[curr_pos].prev_x = x - 1;
        		m[curr_pos].prev_y = y - 1;
        		m[curr_pos].indel_s2 = m[prev_pos].indel_s2;
        		
        		int64_t temp_score;
				/* Tries to use best score from M[x, y-1] in the band*/
        		int64_t left_pos = find_cell(x, y - 1, l1, band);
        		/* If the value is in the band */
				if (left_pos != -1){
        			temp_score = m[left_pos].cell + get_score(s1[y - 1], INDEL);
        			if (m[curr_pos].cell < temp_score){
        				m[curr_pos].cell = temp_score;
        				m[curr_pos].prev_x = x;
        				m[curr_pos].prev_y = y - 1;
        				m[curr_pos].indel_s2 = m[left_pos].indel_s2 + 1;
        			}
        		}
        		
        		/* Tries to use best score from M[x-1, y] in the band*/
        		int64_t up_pos = find_cell(x - 1, y, l1, band);
        		/* If the value is in the band */
        		if (up_pos != -1){
        			temp_score = m[up_pos].cell + get_score(INDEL, s2[x - 1]);
        			if (m[curr_pos].cell < temp_score){
        				m[curr_pos].cell = temp_score;
        				m[curr_pos].prev_x = x - 1;
        				m[curr_pos].prev_y = y;
        				m[curr_pos].indel_s2 = m[up_pos].indel_s2;
					}
				}
			}
		}

    x = l2;
    y = l1;
    /* Start solution's reconstruction */
	int64_t align_length = l2 + m[find_cell(x, y, l1, band)].indel_s2;
	char* s1_align = malloc((align_length + 1) * sizeof(char));
	char* s2_align = malloc((align_length + 1) * sizeof(char));
	char* res_align = malloc((align_length + 1) * sizeof(char));
    s1_align[align_length] = '\0';
    s2_align[align_length] = '\0';
    res_align[align_length] = '\0';
        
	align_s* solution = malloc(1 * sizeof(align_s));
	solution -> score = m[find_cell(x, y, l1, band)].cell;	
		
	int64_t align_index = align_length - 1;
    cell_s c = m[find_cell(x, y, l1, band)];
	/* Reconstruct the solution's strings starting from M[l2, l1] */
	for (; c.prev_x != -1 && c.prev_y != -1;
		x = c.prev_x, y = c.prev_y, c = m[find_cell(x, y, l1, band)]){
		s1_align[align_index] = y == c.prev_y ? '_' : s1[y - 1];
		s2_align[align_index] = x == c.prev_x ? '_' : s2[x - 1];
		if (y == c.prev_y || x == c.prev_x)
			res_align[align_index] = 'I';
		else
			res_align[align_index] = s1[y - 1] == s2[x - 1] ? 'E' : 'U';
		align_index--;
    }
    
	solution -> s1 = s1_align;
	solution -> s2 = s2_align;
	solution -> res = res_align;
	
    return solution;
}

int main(int argc, char **argv) {

	/* Input from stdin */
    int64_t l1, l2;
    
    printf("Insert the length of the first string: \n");
    scanf("%d", &l1);
    assert(l1 > 0 && "length of a string must be > 0");
    char* s1 = malloc((l1 + 1) * sizeof(char));
    printf("Insert the first string: \n");
	scanf("%s", s1);
	       
    printf("Insert the length of the second string: \n");
    scanf("%d", &l2);
    assert(l2 > 0 && "length of a string must be > 0");
    char* s2 = malloc((l2 + 1) * sizeof(char));
    printf("Insert the second string: \n");
	scanf("%s", s2);        
    
	assert(s1 != NULL && "cannot read first string\n");
    assert(s2 != NULL && "cannot read second string\n");
    
    bool swap = false; 
	/* Swap s1 and s2 in order to have s1 as the shortest string and s2 as the longest string */
    if (strlen(s2) < strlen(s1)) {
        char* tmp = s1;
        s1 = s2;
        s2 = tmp;
        swap = true;
    }
	l1 = strlen(s1);
    l2 = strlen(s2);
    align_s* sol_align = NULL;
    /* Band will be the minimum number of indel characters */
	int64_t band = l2 - l1; 
    for (; sol_align == NULL && band / 2 + 1 <= l2 - 1; band *= 2){
        align_s* seq_align = band_align(s1, s2, band);
        int64_t next_bound = maxim_score(band, l2, s1);
#ifdef DEBUG
		printf("Band's size: %d\n", band);
		printf("Best score: %d\n", seq_align -> score);
		printf("Next iteration's upper bound: %d\n", next_bound);
#endif        	
		if (band + 1 > l2 - 1 || next_bound <= seq_align -> score)
			sol_align = seq_align;
	}
    
    /* Output to stdout */
    
    printf("Global alignment score: %d\n", sol_align -> score);
	printf("First string:\n%s\n", swap ? sol_align -> s2 : sol_align -> s1);
	printf("Second string:\n%s\n", swap ? sol_align -> s1 : sol_align -> s2);
	printf("Alignment string (I = Indel, E = Equal, U = Unequal):\n%s\n", sol_align -> res);	
}
