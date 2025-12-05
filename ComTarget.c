
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<ctype.h>
//#include <unistd.h>


#ifdef _WIN32
#include <windows.h>
#include <io.h>  // File Check on Windows
#define F_OK 0
#define access _access
#else
#include <unistd.h>  // File Inspection in Linux/Unix
#include <limits.h>
#endif



/*****************************************************************************************************/
/*                                     ComTarget                                                        */                                                    
/*			Created on Fri December  5 20:46:26 2025                                                      */ 
/*			@author: Yuzhu Li                                                                        */
/*			For more information please contact :	drugyuzhu@163.com                                */
/*													                                                 */                                                                                                     
/*			To compile, use:	gcc ComTarget.c -o ComTarget -lm -O3                                       */
/* 		                                                                                             */						                                                                                                            
/* 		                                                                                             */	          
/* 		                                                                                             */	           
/*****************************************************************************************************/

















struct hu_data {
    char name[100];
    float v1; float v2; float v3; float v4; float v5;
    float v6; float v7; float v8; float v9; float v10;
};

// Comparison Results Structure
struct comparison_result {
    char name[100];
    float similarity;  // Similarity score (closer to 1, more similar)
    int index;         // Index in the original array
};


struct min_max {float xmin;float xmax;float ymin;float ymax;float zmin;float zmax;};
struct grid {float x;float y; float z; int state;short int psa;};
struct data_index{int x;int y;int z;int wsp;};
struct atoms {float xa; float ya; float za;int ele;float rvdw;};
struct atom_info{ int falg; int id; float vdwr; float mw; char atoms[5];};

// Primitive ligand docking scoring structure
struct original_score {
    char name[100];
    float best_score;
    float avg_score;
    float worst_score;
};


// Newly added structure definition
struct score_stats {
    double best_score;    // Best rating (minimum affinity value, strongest combination)
    double average_score; // Average rating
    double worst_score;   // The worst score (the affinity with the highest value, the weakest combination)
    int count;           // Number of modes
};

// Docking result structure
struct docking_result {
    char receptor_name[100];
    struct score_stats docking_scores;  // Current docking statistics information
    struct original_score original_scores;  // Statistical information of the original ligand
    float delta_best;     // Best rating difference
    float delta_avg;      // Average rating difference
    float delta_worst;    // Differential difference
    float delta_composite; // Comprehensive difference (used for sorting)
    int similarity_rank;
};


struct target_info {
    char pdb_id[20];
    char protein_names[500];
    char gene_names[200];
    char uniprot_id[100]; 
};


struct grid *grid_data; 
struct atoms *atom_data;

char* get_valid_filename(int argc, char *argv[]);
int file_exists_and_readable(const char *filename);


struct min_max  get_min_max_mol2(char *name);
struct min_max  get_min_max(struct atoms *atom_data,int atom_number);
struct result{ float surf; float psa; float vol;};
int cal_dis(struct atoms *atom_data,int atom_number,float x,float y,float z);
int cal_adjacent(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number);
int cal_adjacent1(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number);

int cal_psa(struct grid grid_data,struct atoms *atom_data,int atom_number);
int file_type(char *file_name);
int get_atom_info(struct atom_info *atom_info_data);
int get_data(struct atoms **atom_data,struct atom_info *atom_info_data,int count, float x,float y,float z,char *ele_s);
int get_data_from_file(char *file_name,struct atom_info *atom_info_data,struct atoms **atom_data,int *p_atom_number,int *xyz_file_number);
void trace_falef(struct atoms *atom_data_trace,struct atoms *atom_data, int atom_number,int file_number);



int calculate(struct atoms *atom_data, int atom_number,struct hu_data *hudata);
int get_xyz_atom_number(char *name);
int get_xyz_file_number(char *name);
int is_obabel_installed();
int cal_confab(char *file_name);
int cal_conformer(char *file_name);
int cal_obabel_mmf(char *file_name);
int get_mopac_dataf(char *file_name,char *out_name);
int cal_mopac_mmf(char *file_name,int multiplicity,int charge, float permittivity,char *mopac_path );
int get_energy(char *file_name,float bloz);
int get_dist_data(char *file_name,float *distribution, int input_filenumber );
int ask_if_bloz();
float get_dis();
void add_prefix(char *newnamef, char *name, char *add);
void add_suffix(char *newnamef, char *name, char *add);
void set_name(int cal_type,char *input_name,char *seted_name);
void chang_suffix2xyz(char *name);
int is_number(const char *str); 

void cal_3d_hu_moment(struct hu_data *hdata, struct grid  *grid_data, int grid_count);
int cal_moment(struct hu_data *data, struct grid  *grid_data, int grid_count,float grid_gap, float volume, float suface);


int split_xyz(char *file_name);



// Dynamically read descriptor files
int get_discriptorf_dynamic(struct hu_data **descriptors_ptr);

// Calculate the similarity between two descriptors using the division method
float calculate_similarity_ratio(struct hu_data *a, struct hu_data *b);

// Comparison function, used for qsort sorting (descending order, high similarity first)
int compare_similarity_desc(const void *a, const void *b);
// Main function: Compare and sort
int compare_and_sort(struct hu_data *discptor_input, 
                    struct hu_data *discptor_array, 
                    int array_size, 
                    struct comparison_result **results_ptr);



// Release memory function
void free_descriptors(struct hu_data *descriptors) {
    if (descriptors != NULL) {
        free(descriptors);
    }
}

void free_results(struct comparison_result *results) {
    if (results != NULL) {
        free(results);
    }
}



int gen_conf(char *input_file) ;


void compute_descriptor(char *filename, struct atom_info *atom_info_data, struct hu_data *output_descriptor);



int covert2pdbqt(char *name);



int get_boxf(char *name, char *box_data);
//void rundockf(char *receptor,char *ligand,int id);
struct score_stats rundockf(char *receptor, char *ligand, int id);

int dock_infof(int max_targets);

// Comparison function sorted by comprehensive difference (ascending order, the smaller the difference, the better)
int compare_delta_composite_asc(const void *a, const void *b); 
// Comparison function sorted by the best docking score (ascending order, the smaller the score, the better)
int compare_docking_best_asc(const void *a, const void *b);

// Read the original ligand docking scoring file
int read_original_scores(struct original_score **scores_ptr);

int extract_scores(const char *filename, struct score_stats *stats);

int preview_boltzmann_distribution(char *file_name); 
float get_dis_with_preview(char *file_name);
int compare_delta_best_asc(const void *a, const void *b);


//Read the Targets. tsv file 
//Declaration of target information related functions 
//Extract PDB ID from composite name
void extract_pdb_id(const char *compound_name, char *pdb_id); 
int read_targets_info(const char *filename, struct target_info **targets_ptr);
void extract_pdb_id(const char *compound_name, char *pdb_id);
struct target_info* get_target_info(const char *compound_name, struct target_info *targets, int target_count);
void print_single_target_info(int rank, const char *pdb_id, struct target_info *target);
char* string(int length, char character);
void free_targets_info(struct target_info *targets);



void save_results_to_csv(struct docking_result *docking_results, int dock_num, 
                        struct target_info *targets, int target_count, 
                        const char *filename); 


void print_smart_progress(int current, int total); 

char* get_current_time();


void print_progress_bar1(int current, int total, const char* description); 

int get_rotatable_bonds(char *file_name);

void prompt_delete_temp_files();
int compare_docking_best_then_delta(const void *a, const void *b);

void get_program_dir(char *buffer, size_t size);


int main(int argc, char *argv[])
{

 
	struct hu_data input_descriptor;	
	int atom_number,i,file_number,if_trace_file=0;
	struct atom_info *atom_info_data;
	atom_info_data=(struct atom_info *)malloc(200*sizeof(struct atom_info));
	get_atom_info(atom_info_data);


	/**
		**Read info from the input file
	**/	
	char *file_name=NULL;
	//file_name=argv[1];
	file_name = get_valid_filename(argc, argv);
	if (file_name == NULL) {
        printf("Unable to obtain valid files, program exits\n");
        free(atom_info_data);
        return -1;
    }
	
	
	
	
	int filenumber;
	
	filenumber =gen_conf(file_name);
	
    if (filenumber == 0) {
        printf("Conformation generation failed\n");
        free(atom_info_data);
        return -1;
    }	
	printf("Successfully generated %d conformational files\n", filenumber);
	
    struct hu_data *descriptors = NULL;
    int count = get_discriptorf_dynamic(&descriptors);	
    printf("Successfully read %d database descriptors\n", count);
    
    if (count == 0) {
        printf("No descriptor data was read\n");
        free(atom_info_data);
        return -1;
    }
	
    // Create an array to store all results
    int max_total_results = filenumber * 200; // Each file can retrieve up to 100 results
    struct comparison_result *all_results = (struct comparison_result *)malloc(max_total_results * sizeof(struct comparison_result));
    int total_results_count = 0;
	
	if(filenumber>1)
	{ 
		// Process each generated conformation file
   	 	for (i = 1; i <= filenumber; i++) {
     	   char temp_filename[50];
      	  sprintf(temp_filename, "temp_%d.xyz", i);
        
      	  printf("Process files: %s\n", temp_filename);
        
      	  // Calculate the descriptor of the current file
      	  struct hu_data current_descriptor;
      	  compute_descriptor(temp_filename, atom_info_data, &current_descriptor);
        
      	  // Compare with database
      	  struct comparison_result *current_results = NULL;
      	  int current_count = compare_and_sort(&current_descriptor, descriptors, count, &current_results);
        
      	  if (current_count == 0) {
          	  printf("File %s comparison failed\n", temp_filename);
          	  continue;
      	  }
        
      	  // Add the top 100 results to the summary set
      	  int take_count = (current_count < 200) ? current_count : 200;
      	  int j;
       	 for ( j = 0; j < take_count; j++) {
           	 if (total_results_count < max_total_results) {
           	     all_results[total_results_count] = current_results[j];
          	      total_results_count++;
          	  }
      	  }
        
       	 free_results(current_results);
   	 	}		
	}	
	if(filenumber==1)
	{
		// Process each generated conformation file
   	 	for (i = 1; i <= filenumber; i++) {
   	 		
     		char temp_filename[50],command[300];
     		command[0]='\0';
     		strcat(command,"obabel ");
     		strcat(command,file_name);
     		strcat(command," -o xyz -O ");
     		sprintf(temp_filename, "%s",file_name);
     		int j=0;
     		while(temp_filename[j]!='\0') j++;
     		while(temp_filename[j]!='.') j--;
     		temp_filename[j]='\0';
     		strcat(temp_filename,".xyz");
     		strcat(command,temp_filename);
     		strcat(command,"  -d");
     		system(command);
      	  printf("Process files: %s\n", temp_filename);
        
      	  // Calculate the descriptor of the current file
      	  struct hu_data current_descriptor;
      	  compute_descriptor(temp_filename, atom_info_data, &current_descriptor);
        
      	  // Compare with database
      	  struct comparison_result *current_results = NULL;
      	  int current_count = compare_and_sort(&current_descriptor, descriptors, count, &current_results);
        
      	  if (current_count == 0) {
          	  printf("File %s comparison failed\n", temp_filename);
          	  continue;
      	  }
        
      	  // Add the top 200 results to the summary set
      	  int take_count = (current_count < 200) ? current_count : 200;
      	 
       	 for ( j = 0; j < take_count; j++) {
           	 if (total_results_count < max_total_results) {
           	     all_results[total_results_count] = current_results[j];
          	      total_results_count++;
          	  }
      	  }
        
       	 free_results(current_results);
   	 	}		
	}

    // If there is no result, exit
    if (total_results_count == 0) {
        printf("No similar results were found\n");
        free(all_results);
        free_descriptors(descriptors);
        free(atom_info_data);
        return -1;
    }

    // Sort all results in descending order of similarity
    qsort(all_results, total_results_count, sizeof(struct comparison_result), compare_similarity_desc);

    // De duplication: Retain the highest similarity result of each molecule
    struct comparison_result *unique_results = (struct comparison_result *)malloc(total_results_count * sizeof(struct comparison_result));
    int unique_count = 0;
    
    // Used to record the names of molecules that have appeared
    char **seen_names = (char **)malloc(total_results_count * sizeof(char *));
    int seen_count = 0;
    
    for (i = 0; i < total_results_count; i++) {
        int found = 0,j;
        
        // Check if this molecule has been processed before
        for ( j = 0; j < seen_count; j++) {
            if (strcmp(all_results[i].name, seen_names[j]) == 0) {
                found = 1;
                break;
            }
        }
        
        // If it is a new molecule, add it to the unique result list
        if (!found) {
            unique_results[unique_count] = all_results[i];
            seen_names[seen_count] = strdup(all_results[i].name);
            seen_count++;
            unique_count++;
        }
    }

	int actual_unique_count = unique_count;

    // Create similarity result file name
    char similarity_filename[256];
    char temp_name[256];
    sprintf(temp_name, "%s", file_name); 
    i=0;
    while(temp_name[i]!='\0') i++;
    while(temp_name[i]!='.') i--;
    temp_name[i]='\0';
   // sprintf(similarity_filename, "similarity_results_%s.txt", temp_name);
    snprintf(similarity_filename, sizeof(similarity_filename), "similarity_results_%s.txt", temp_name);
    // Open the similarity result file
    FILE *similarity_file = fopen(similarity_filename, "w");
    if (similarity_file == NULL) {
        printf("Warning: Could not open %s for writing. Similarity results will not be saved.\n", similarity_filename);
    } else {
        fprintf(similarity_file, "Similarity Analysis Results - Generated on %s\n", get_current_time());
        fprintf(similarity_file, "Input file: %s\n\n", file_name);
        fflush(similarity_file);
    }

	// Output the top 200 most similar results

    if (similarity_file) {
        fprintf(similarity_file, "\nFinal result (top %d most similar descriptors after deduplication):\n", 
                (unique_count < 200) ? unique_count : 200);
        int print_count = (unique_count < 200) ? unique_count : 200;
        for (i = 0; i < print_count; i++) {
            fprintf(similarity_file, "%d. %s (similarity: %.4f)\n", 
                   i + 1, unique_results[i].name, unique_results[i].similarity);
        }
        fprintf(similarity_file, "\nTotal unique results found: %d\n", unique_count);
        fclose(similarity_file);
        printf("Similarity analysis completed! Results saved to: %s\n", similarity_filename);
    }


    
 /*   printf("\nFinal result (top 100 most similar descriptors after deduplication):\n");
    int print_count = (unique_count < 200) ? unique_count : 200;
    for (i = 0; i < print_count; i++) {
        printf("%d. %s (similarity: %.4f)\n", 
               i + 1, unique_results[i].name, unique_results[i].similarity);
    }  */

/***************************************************************/
	// Docking calculation
	
	
covert2pdbqt(file_name);
int dock_num;
dock_num = dock_infof(actual_unique_count);

if (dock_num > 0) {
    // Read the original ligand score
    struct original_score *original_scores = NULL;
    int original_count = read_original_scores(&original_scores);
    
    if (original_count == 0) {
        printf("Warning: No original score data found, will use docking scores only\n");
    }
    
    // Create docking result array
    struct docking_result *docking_results = (struct docking_result *)malloc(dock_num * sizeof(struct docking_result));
    
   // printf("\nPerforming docking calculations...\n");
   
   
	printf("\n\n"
		" _______________________________________________________________________ \n"
		"|                                                                       |\n"
		"|                        Reverse Docking Section                        |\n"
		"|_______________________________________________________________________|\n\n");	
	   
   
    
        printf("\nPerforming docking calculations...\n");
        printf("Progress: 0%%   10   20   30   40   50   60   70   80   90   100%%\n");
        printf("         |----|----|----|----|----|----|----|----|----|----|\n");
        printf("         ");
        fflush(stdout);
        // Initialization progress tracking
        int last_percentage = -1;        
    
// Open the file to save detailed information
char detail_filename[256];
snprintf(detail_filename, sizeof(detail_filename), "docking_details_%s.txt", temp_name);
FILE *detail_file = fopen(detail_filename, "w");

if (detail_file == NULL) {
    printf("Warning: Could not open %s for writing. Details will not be saved.\n", detail_filename);
} else {
    // Add timestamp and basic information at the beginning of the file
    fprintf(detail_file, "Docking Details - Generated on %s\n", get_current_time());
    fprintf(detail_file, "Input file: %s\n", file_name);
    fprintf(detail_file, "Number of targets: %d\n\n", dock_num);
    fflush(detail_file);
}

for (i = 0; i < dock_num; i++) {
    strcpy(docking_results[i].receptor_name, unique_results[i].name);
    docking_results[i].similarity_rank = i + 1;
    
    // Run docking and obtain complete scoring statistics
    docking_results[i].docking_scores = rundockf(unique_results[i].name, file_name, i + 1);
    
    // Find the corresponding original ligand score
    int found = 0;
    int j;
    for (j = 0; j < original_count; j++) {
        if (strcmp(original_scores[j].name, unique_results[i].name) == 0) {
            docking_results[i].original_scores = original_scores[j];
            found = 1;
            break;
        }
    }
    
    if (!found) {
        printf("Warning: No original score found for %s\n", unique_results[i].name);
        // Set default values
        docking_results[i].original_scores.best_score = 0.0;
        docking_results[i].original_scores.avg_score = 0.0;
        docking_results[i].original_scores.worst_score = 0.0;
    }
    
    // Calculate various differences and handle missing data (1000 indicates missing)
    if (docking_results[i].original_scores.best_score >= 1000.0f) {
        // The original ligand data is missing, and the difference will not be calculated
        docking_results[i].delta_best = 999.0f; // Special markings indicate missing information
    } else {
        docking_results[i].delta_best = docking_results[i].docking_scores.best_score - 
                                       docking_results[i].original_scores.best_score;
    }
    
    if (docking_results[i].original_scores.avg_score >= 1000.0f) {
        docking_results[i].delta_avg = 999.0f;
    } else {
        docking_results[i].delta_avg = docking_results[i].docking_scores.average_score - 
                                      docking_results[i].original_scores.avg_score;
    }
    
    if (docking_results[i].original_scores.worst_score >= 1000.0f) {
        docking_results[i].delta_worst = 999.0f;
    } else {
        docking_results[i].delta_worst = docking_results[i].docking_scores.worst_score - 
                                        docking_results[i].original_scores.worst_score;
    }
    
    // Calculate the comprehensive difference and handle missing data
    if (docking_results[i].delta_best >= 999.0f || 
        docking_results[i].delta_avg >= 999.0f || 
        docking_results[i].delta_worst >= 999.0f) {
        // If any data is missing, the comprehensive difference is also marked as missing
        docking_results[i].delta_composite = 999.0f;
    } else {
        // Calculate the comprehensive difference (weight: best 50%, average 30%, worst 20%)
        docking_results[i].delta_composite = (docking_results[i].delta_best * 0.5f) +
                                           (docking_results[i].delta_avg * 0.3f) +
                                           (docking_results[i].delta_worst * 0.2f);
    }

    // Write detailed information to a file (not displayed on the terminal)
    if (detail_file != NULL) {
        fprintf(detail_file, "Docking %d/%d: %s\n", i + 1, dock_num, unique_results[i].name);
        fprintf(detail_file, "   Current - Best: %.2f, Avg: %.2f, Worst: %.2f\n",
               docking_results[i].docking_scores.best_score,
               docking_results[i].docking_scores.average_score,
               docking_results[i].docking_scores.worst_score);
        
        // Check if the original data is missing
        if (docking_results[i].original_scores.best_score >= 1000.0f) {
            fprintf(detail_file, "   Original - [DATA MISSING]\n");
        } else {
            fprintf(detail_file, "   Original - Best: %.2f, Avg: %.2f, Worst: %.2f\n",
                   docking_results[i].original_scores.best_score,
                   docking_results[i].original_scores.avg_score,
                   docking_results[i].original_scores.worst_score);
        }
        
        // Check if the difference is missing
        if (docking_results[i].delta_best >= 999.0f) {
            fprintf(detail_file, "   Delta - [DATA MISSING]\n");
        } else {
            fprintf(detail_file, "   Delta - Best: %.2f, Avg: %.2f, Worst: %.2f, Composite: %.2f\n",
                   docking_results[i].delta_best, docking_results[i].delta_avg,
                   docking_results[i].delta_worst, docking_results[i].delta_composite);
        }
        fprintf(detail_file, "\n");
        fflush(detail_file);
    }
    
    // Only display progress bar, do not show detailed information
    print_smart_progress(i+1, dock_num);
}

// Close detailed files
if (detail_file != NULL) {
    fclose(detail_file);
    printf("\n\n Docking completed! All %d targets processed.\n", dock_num);
    printf("   Detailed results saved to: %s\n", detail_filename);
} else {
    printf("\n\n Docking completed! All %d targets processed.\n", dock_num);
}

// Read target information
struct target_info *targets = NULL;
int target_count = read_targets_info("./info_data/Targets.tsv", &targets);

if (target_count == 0) {
    printf("Warning: Could not load target information from Targets.tsv\n");
    printf("Trying alternative path: Targets.tsv\n");
    target_count = read_targets_info("Targets.tsv", &targets);
}

// Use a new sorting function: First, sort by the best docking score, and for those with the same score, sort by the difference in the best score
//printf("\n=== Final Prediction Results (Primary Sort by Best Docking Score) ===\n");
qsort(docking_results, dock_num, sizeof(struct docking_result), compare_docking_best_then_delta);
/*
for (i = 0; i < dock_num; i++) {
    printf("%d. Target: %s\n", i + 1, docking_results[i].receptor_name);
    printf("   Best Score: %.2f kcal/mol", docking_results[i].docking_scores.best_score);
    
    // Display difference information (if available)
    if (docking_results[i].delta_best < 999.0f) {
        printf(" (Delta: %.2f)", docking_results[i].delta_best);
    } else {
        printf(" [Original data missing]");
    }
    printf("\n");
    
    // Simplified prediction suggestions for handling missing data
    if (docking_results[i].delta_best >= 999.0f) {
        printf("   Prediction: DATA MISSING - cannot evaluate\n");
    } else if (docking_results[i].delta_best < -1.0) {
        printf("   Prediction: HIGH potential target\n");
    } else if (docking_results[i].delta_best < -0.5) {
        printf("   Prediction: GOOD potential target\n");
    } else if (docking_results[i].delta_best < 0.5) {
        printf("   Prediction: MEDIUM potential target\n");
    } else if (docking_results[i].delta_best < 1.5) {
        printf("   Prediction: LOW potential target\n");
    } else {
        printf("   Prediction: VERY LOW potential target\n");
    }
    printf("   --------------------\n");
}   */

// Display target information
//printf("\n=== Target Information (Sorted by Prediction Rank) ===\n");
//printf("This section provides detailed biological information for each predicted target.\n");
//printf("%s\n", string(80, '='));

for (i = 0; i < dock_num; i++) {
    struct target_info *target = get_target_info(docking_results[i].receptor_name, targets, target_count);
 //   print_single_target_info(i + 1, docking_results[i].receptor_name, target);
}

// Save results to CSV file
char docking_results_name[200];
snprintf(docking_results_name, sizeof(docking_results_name), "docking_results_%s.csv", temp_name);
save_results_to_csv(docking_results, dock_num, targets, target_count, docking_results_name);

prompt_delete_temp_files();
 
    
    // release memory
    free_targets_info(targets);
    free(docking_results);
    if (original_scores != NULL) {
        free(original_scores);
    }
}	
	
	
	
	

	
	
	
	
	
	
	
	
	

    // release memory
    free(all_results);
    for (i = 0; i < seen_count; i++) {
        free(seen_names[i]);
    }
    free(seen_names);
    free(unique_results);
    free_descriptors(descriptors);
    free(atom_info_data);

   


return 0;
}




void get_program_dir(char *buffer, size_t size) {
    #ifdef _WIN32
        char exe_path[MAX_PATH];
        GetModuleFileName(NULL, exe_path, MAX_PATH);
        // 找到最后一个路径分隔符
        char *last_slash = strrchr(exe_path, '\\');
        char *last_slash2 = strrchr(exe_path, '/');
        if (last_slash2 > last_slash) last_slash = last_slash2;
        if (last_slash) {
            size_t len = last_slash - exe_path;
            if (len < size) {
                strncpy(buffer, exe_path, len);
                buffer[len] = '\0';
            } else {
                buffer[0] = '\0';
            }
        } else {
            buffer[0] = '\0';
        }
    #else
        char exe_path[PATH_MAX];
        ssize_t count = readlink("/proc/self/exe", exe_path, sizeof(exe_path)-1);
        if (count != -1) {
            exe_path[count] = '\0';
            char *last_slash = strrchr(exe_path, '/');
            if (last_slash) {
                size_t len = last_slash - exe_path;
                if (len < size) {
                    strncpy(buffer, exe_path, len);
                    buffer[len] = '\0';
                } else {
                    buffer[0] = '\0';
                }
            } else {
                buffer[0] = '\0';
            }
        } else {
            buffer[0] = '\0';
        }
    #endif
}
















// New comparison function: First, sort by the best docking score, and for those with the same score, sort by the difference in the best score

int compare_docking_best_then_delta(const void *a, const void *b) {
    struct docking_result *resultA = (struct docking_result *)a;
    struct docking_result *resultB = (struct docking_result *)b;
    
    // Firstly, compare the best docking scores (in ascending order, the smaller the score, the better)
    if (resultA->docking_scores.best_score < resultB->docking_scores.best_score) return -1;
    if (resultA->docking_scores.best_score > resultB->docking_scores.best_score) return 1;
    
    // If the best scores are the same, compare the difference between the best scores (in ascending order, the smaller the difference, the better)
    if (resultA->delta_best < resultB->delta_best) return -1;
    if (resultA->delta_best > resultB->delta_best) return 1;
    
    return 0;
}

/**
 * Check if the file exists and is readable
 * Return: 1- Exists and is readable, 0- Does not exist or is unreadable
 */
int file_exists_and_readable(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (file != NULL) {
        fclose(file);
        return 1;
    }
    return 0;
}

/**
 * Get a valid file name
 * Return: A valid file name pointer, returns NULL if failed
 */
char* get_valid_filename(int argc, char *argv[])
{
    static char filename[256];
    int max_attempts = 3;
    int attempts = 0;
    
    while (attempts < max_attempts) {
        // Check command line parameters
        if (argc < 2 || attempts > 0) {
            if (attempts == 0) {
                printf("Please enter the file path: ");
            } else {
                printf("Please re-enter the file path: ");
            }
            
            if (fgets(filename, sizeof(filename), stdin) == NULL) {
                printf("Read input failed\n");
                return NULL;
            }
            
            // Remove line breaks
            filename[strcspn(filename, "\n")] = 0;
        } else {
            strncpy(filename, argv[1], sizeof(filename)-1);
            filename[sizeof(filename)-1] = '\0';
        }
        
        // Check if the file is an empty string
        if (strlen(filename) == 0) {
            printf("Error: The file path cannot be empty\n");
            attempts++;
            continue;
        }
        
        // Check if the file exists and is readable
        if (!file_exists_and_readable(filename)) {
            printf("Error: file '%s' Does not exist or is unreadable\n", filename);
            attempts++;
            continue;
        }
        
        // Valid document
        return filename;
    }
    
    printf("Exceeding the maximum number of retries (%d)，Program exit\n", max_attempts);
    return NULL;
}







void prompt_delete_temp_files() {
    int opt = 0;
    char buffer_command[300];
    int valid_input = 0; // Mark whether the input is valid

    printf("Congratulations, all calculations have been completed!\n");
    
    while (!valid_input) {
        printf("Do you need to delete temporary files generated during the calculation process? 1 (Yes) / 0 (No)\n");
        printf("Please enter your choice: ");
        
        if (fgets(buffer_command, sizeof(buffer_command), stdin) == NULL) {
            printf("Error reading input. Exiting.\n");
            return;
        }
        
        // Check if the input is an integer
        if (sscanf(buffer_command, "%d", &opt) == 1) {
            if (opt == 1) {
                #ifdef _WIN32
                int ret1 = system("del temp*.xyz");
                int ret2 = system("del *.log");
                #else
                int ret1 = system("rm temp*.xyz");
                int ret2 = system("rm *.log");
                #endif
                // Optional: Check command execution results
                if (ret1 != 0 || ret2 != 0) {
                    printf("Warning: Some files may not have been deleted successfully.\n");
                } else {
                    printf("Temporary files temp*.xyz and *.log have been deleted.\n");
                }
                valid_input = 1;
            } else if (opt == 0) {
                printf("Temporary files preserved.\n");
                valid_input = 1;
            } else {
                printf("Invalid input. Please enter 1 (Yes) or 0 (No).\n");
            }
        } else {
            printf("Invalid input. Please enter a number (1 or 0).\n");
        }
    }
}











int get_rotatable_bonds(char *file_name) {
    char command[300];
    char output_filename[100];
    FILE *fp;
    int rot_bonds = -1;

    // Create a unique temporary file name
    sprintf(output_filename, "temp_rotors_%d.txt", (int)time(NULL));
    
    // Build and execute commands
    sprintf(command, "obabel %s -o smi --append \"rotors\" > %s 2>&1", file_name, output_filename);
    system(command);
 
    fp = fopen(output_filename, "r");
    if (fp != NULL) {
        char output[1000];
        int line_count = 0;
        
        // Read all rows
        while (fgets(output, sizeof(output), fp) != NULL) {
            line_count++;
            // The second line contains molecular information (skipping the state information in the first line)
            if (line_count == 2) {              
                // Find the last space
                char *last_space = strrchr(output, ' ');
                if (last_space != NULL) {
                    char *rotor_str = last_space + 1;
                    rotor_str[strcspn(rotor_str, "\n\r")] = '\0';
                    rot_bonds = atoi(rotor_str);
                    break;
                }
            }
        }
        fclose(fp);
        
        if (line_count < 2) {
            printf("Error: Expected at least 2 lines in output, but got %d\n", line_count);
        }
    }
    
    // Clean up temporary files
    remove(output_filename);
    
    return rot_bonds;
}




















void print_progress_bar1(int current, int total, const char* description) {
    int barWidth = 50;
    float progress = (float)current / total;
    int pos = barWidth * progress;
    int i;

    printf("\r%s [", description);
    for (i = 0; i < barWidth; ++i) {
        if (i < pos) printf("=");
        else if (i == pos) printf(">");
        else printf(" ");
    }
    printf("] %d%% (%d/%d)", (int)(progress * 100), current, total);
    fflush(stdout);
}












char* get_current_time() {
    static char buffer[80];
    time_t rawtime;
    struct tm *timeinfo;
    
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);
    return buffer;
}










void print_smart_progress(int current, int total) {
    const int total_width = 51; // 进度条总宽度，与框架字符数一致
    int stars; // 当前应显示的星号数量
    static int last_stars = 0,i; // 上次打印的星号数量

    // 计算当前星号数量，使用向上取整以避免进度滞后
    stars = (current * total_width + total - 1) / total;

    // 确保 stars 在合理范围内
    if (stars < 0) stars = 0;
    if (stars > total_width) stars = total_width;

    // 打印缺失的星号（从上次位置到当前位置）
    if (stars > last_stars) {
        for ( i = last_stars; i < stars; i++) {
            printf("*");
        }
        fflush(stdout); // 立即输出，避免缓冲
        last_stars = stars;
    }

    // 最终检查：当完成时，确保打印所有星号（防止舍入误差）
    if (current == total && last_stars < total_width) {
        for ( i = last_stars; i < total_width; i++) {
            printf("*");
        }
        fflush(stdout);
        last_stars = total_width;
    }
}





// 按最佳打分差值排序的比较函数（升序，差值越小越好）
int compare_delta_best_asc(const void *a, const void *b) {
    struct docking_result *resultA = (struct docking_result *)a;
    struct docking_result *resultB = (struct docking_result *)b;
    
    if (resultA->delta_best < resultB->delta_best) return -1;
    if (resultA->delta_best > resultB->delta_best) return 1;
    return 0;
}






int extract_scores(const char *filename, struct score_stats *stats) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }

    //Initialize statistics 
    int DBL_MAX = 1000;
    stats->best_score = DBL_MAX;    //Find the minimum value (the most negative value) 
    stats->worst_score = -DBL_MAX;  //Find the maximum value (the least negative value) 
    stats->average_score = 0.0;
    stats->count = 0;
    
    double sum_affinity = 0.0;
    char line[256];

    //Read files line by line 
    while (fgets(line, sizeof(line), file)) {
        int mode;
        double affinity, rmsd_lb, rmsd_ub;
        
        //Attempt to parse data rows 
        if (sscanf(line, "%d %lf %lf %lf", &mode, &affinity, &rmsd_lb, &rmsd_ub) == 4) {
            //Update the best rating (search for smaller values, i.e. more negative values) 
            if (affinity < stats->best_score) {
                stats->best_score = affinity;
            }
            //Update the worst score (search for a larger value, i.e. a more non negative value) 
            if (affinity > stats->worst_score) {
                stats->worst_score = affinity;
            }
            sum_affinity += affinity;
            stats->count++;
        }
    }

    fclose(file);

    //Check if data has been found 
    if (stats->count == 0) {
        fprintf(stderr, "No data lines found in file\n");
        return -1;
    }

    //Calculate the average affinity
    stats->average_score = sum_affinity / stats->count;

    return 0;
}









 











// Comparison function, used for qsort sorting (descending order, high similarity first)
int compare_similarity_desc(const void *a, const void *b) {
    struct comparison_result *resultA = (struct comparison_result *)a;
    struct comparison_result *resultB = (struct comparison_result *)b;
    
    if (resultA->similarity > resultB->similarity) return -1;
    if (resultA->similarity < resultB->similarity) return 1;
    return 0;
}

// Main function: Compare and sort
int compare_and_sort(struct hu_data *discptor_input, 
                    struct hu_data *discptor_array, 
                    int array_size, 
                    struct comparison_result **results_ptr) {
    
    if (discptor_input == NULL || discptor_array == NULL || results_ptr == NULL) {
        return 0;
    }
    
    // Allocation result array
    *results_ptr = (struct comparison_result*)malloc(array_size * sizeof(struct comparison_result));
    if (*results_ptr == NULL) {
        printf("Result array memory allocation failed\n");
        return 0;
    }
    
    // Calculate the similarity between each descriptor and the input descriptor
    int i;
    for ( i = 0; i < array_size; i++) {
        strcpy((*results_ptr)[i].name, discptor_array[i].name);
        (*results_ptr)[i].similarity = calculate_similarity_ratio(discptor_input, &discptor_array[i]);
        (*results_ptr)[i].index = i;
    }
    
    // Sort in descending order of similarity (the larger the similarity value, the more similar it is)
    qsort(*results_ptr, array_size, sizeof(struct comparison_result), compare_similarity_desc);
    
    return array_size;
}









// Calculate the similarity between two descriptors using the division method
float calculate_similarity_ratio(struct hu_data *a, struct hu_data *b) {
    float similarity = 0.0;
    int valid_dimensions = 0;
    
    // Calculate the ratio for each dimension and add it to the similarity score
    if (b->v1 != 0) {
        float ratio = a->v1 / b->v1;
        // The closer the ratio is to 1, the more similar it is. Use 1- | 1-ratio | to calculate similarity
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v2 != 0) {
        float ratio = a->v2 / b->v2;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v3 != 0) {
        float ratio = a->v3 / b->v3;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v4 != 0) {
        float ratio = a->v4 / b->v4;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v5 != 0) {
        float ratio = a->v5 / b->v5;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v6 != 0) {
        float ratio = a->v6 / b->v6;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v7 != 0) {
        float ratio = a->v7 / b->v7;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v8 != 0) {
        float ratio = a->v8 / b->v8;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v9 != 0) {
        float ratio = a->v9 / b->v9;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    if (b->v10 != 0) {
        float ratio = a->v10 / b->v10;
        similarity += 1.0 - fabsf(1.0 - ratio);
        valid_dimensions++;
    }
    
    // Return the average similarity
    return (valid_dimensions > 0) ? (similarity / valid_dimensions) : 0.0;
}



// Dynamically read descriptor files
int get_discriptorf_dynamic(struct hu_data **descriptors_ptr) {
    FILE *file = fopen("./info_data/discriptor.txt", "r");
    if (file == NULL) {
		printf("Error: Could not open descriptor file './info_data/discriptor.txt'\n");
		printf("Please ensure the file exists in the correct directory.\n");
		perror("File open error");        
        *descriptors_ptr = NULL;
        return 0;
    }
    
    // First time: Calculate the number of rows
    int count = 0;
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        count++;
    }
    
    // allocate memory
    *descriptors_ptr = (struct hu_data*)malloc(count * sizeof(struct hu_data));
    if (*descriptors_ptr == NULL) {
        printf("allocate memory failed\n");
        fclose(file);
        return 0;
    }
    
    // Second pass: Read data
    rewind(file);
    int actual_count = 0;
    
    while (fgets(line, sizeof(line), file)) {
        // Remove line breaks at the end of a line
        line[strcspn(line, "\n")] = 0;
        
        // Skip blank lines
        if (strlen(line) == 0) continue;
        
        // Using strtok to split strings
        char *token = strtok(line, "#");
        if (token == NULL) continue;
        
        // Copy Name
        strncpy((*descriptors_ptr)[actual_count].name, token, 
                sizeof((*descriptors_ptr)[actual_count].name) - 1);
        (*descriptors_ptr)[actual_count].name[sizeof((*descriptors_ptr)[actual_count].name) - 1] = '\0';
        
        // Read 10 floating-point values
        int value_count = 0;
        float values[10];
        
        while ((token = strtok(NULL, "#")) != NULL && value_count < 10) {
            values[value_count] = atof(token);
            value_count++;
        }
        
        // If 10 values are successfully read, fill the structure
        if (value_count == 10) {
            (*descriptors_ptr)[actual_count].v1 = values[0];
            (*descriptors_ptr)[actual_count].v2 = values[1];
            (*descriptors_ptr)[actual_count].v3 = values[2];
            (*descriptors_ptr)[actual_count].v4 = values[3];
            (*descriptors_ptr)[actual_count].v5 = values[4];
            (*descriptors_ptr)[actual_count].v6 = values[5];
            (*descriptors_ptr)[actual_count].v7 = values[6];
            (*descriptors_ptr)[actual_count].v8 = values[7];
            (*descriptors_ptr)[actual_count].v9 = values[8];
            (*descriptors_ptr)[actual_count].v10 = values[9];
            
            actual_count++;
        }
    }
    
    fclose(file);
    return actual_count;
}











//Comparison function sorted by comprehensive difference (ascending order, the smaller the difference, the better)
int compare_delta_composite_asc(const void *a, const void *b) {
    struct docking_result *resultA = (struct docking_result *)a;
    struct docking_result *resultB = (struct docking_result *)b;
    
    if (resultA->delta_composite < resultB->delta_composite) return -1;
    if (resultA->delta_composite > resultB->delta_composite) return 1;
    return 0;
}

//Comparison function sorted by the best docking score (ascending order, the smaller the score, the better)
int compare_docking_best_asc(const void *a, const void *b) {
    struct docking_result *resultA = (struct docking_result *)a;
    struct docking_result *resultB = (struct docking_result *)b;
    
    if (resultA->docking_scores.best_score < resultB->docking_scores.best_score) return -1;
    if (resultA->docking_scores.best_score > resultB->docking_scores.best_score) return 1;
    return 0;
}













//Read the original ligand docking scoring file
int read_original_scores(struct original_score **scores_ptr) {
    FILE *file = fopen("./info_data/score_data.txt", "r");
    if (file == NULL) {
        printf("Could not open file score_data.txt\n");
        *scores_ptr = NULL;
        return 0;
    }
    
    //Calculate the number of rows for the first time
    int count = 0;
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        count++;
    }
    
    // allocate memory
    *scores_ptr = (struct original_score*)malloc(count * sizeof(struct original_score));
    if (*scores_ptr == NULL) {
        printf("Memory allocation failed for original scores\n");
        fclose(file);
        return 0;
    }
    
    // Second reading of data
    rewind(file);
    int actual_count = 0;


 



    
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;
        
        //Skip blank lines 
        if (strlen(line) == 0) continue;
        
        //Using strtok to split strings
        char *token = strtok(line, "#");
        if (token == NULL) continue;
        
        //Copy Name 
        strncpy((*scores_ptr)[actual_count].name, token, 
                sizeof((*scores_ptr)[actual_count].name) - 1);
        (*scores_ptr)[actual_count].name[sizeof((*scores_ptr)[actual_count].name) - 1] = '\0';
        
        //Read 3 rating points
        token = strtok(NULL, "#");
        if (token == NULL) continue;
        (*scores_ptr)[actual_count].best_score = atof(token);
        
        token = strtok(NULL, "#");
        if (token == NULL) continue;
        (*scores_ptr)[actual_count].avg_score = atof(token);
        
        token = strtok(NULL, "#");
        if (token == NULL) continue;
        (*scores_ptr)[actual_count].worst_score = atof(token);
        
        actual_count++;
    }
    
    fclose(file);
    return actual_count;
}



















	int dock_infof(int max_targets)
	{
		
		char buffer_command[300],newname[100],newname1[100];
		int dock_sel=0,dock_num=0;
		printf("\n Do you need to use molecular docking for secondary screening of similarity calculation results? \n\nYes(1), not(0):\n");		
		if(fgets(buffer_command, sizeof(buffer_command), stdin)==NULL){
			fprintf(stderr,"Error reading input.\n");
			return 1;
		}
		sscanf(buffer_command, "%d", &dock_sel);	
		if(dock_sel==1)
		{
			printf("You chooed method is：[1] Doncking will by Autodock Vina\n\n");
			printf("Please select the number of docking required (1-%d)\n", max_targets);
			
			if(fgets(buffer_command, sizeof(buffer_command), stdin)==NULL){
				fprintf(stderr,"Error reading input.\n");
				return 0;
			}
			sscanf(buffer_command, "%d", &dock_num);			
			if((dock_num>=1)&&(dock_num<=max_targets)){
				printf("\n\nYou will perform docking calculations based on the top %d similarity filters.\n",dock_num);
				return dock_num;
				
			} 
			else{
				fprintf(stderr,"Error reading input.\n");
				return 0;
			}
			
			
			
		} 	
		else if(dock_sel == 0) return 0;		
		
		
		
		
		
	}




int get_boxf(char *name, char *box_data)
{
	char temp_name[100];
	temp_name[0]='#';
	temp_name[1]='\0';
	strcat(temp_name,name);
	FILE *boxfp;
	int i,strsize;
	char fstring[200];
	boxfp=fopen("./info_data/box.data","r");
	strsize=strlen(temp_name);
	if(boxfp!=NULL)
	{

		while((fgets(fstring,199,boxfp))!=NULL)
		{
			if((strncmp(fstring,temp_name,strsize)==0)&&(fstring[strsize]=='#'))
			{
				//printf("%s\n",fstring);
				box_data[0]='\0';
				strcat(box_data," --receptor ./info_data/receptor_data/");
				int j=0;
				i=0;
				char temmp[5];
				for(i=3;i<=6;i++)
				{
					temmp[j]=fstring[i];
				
					j++;
				 } 
				temmp[j]='\0';
				strcat(box_data,temmp);
				strcat(box_data,".pdbqt  ");
				i=0;
				fstring[i]=' ';
				while(fstring[i]!='#')
				{
					fstring[i]=' ';
					i++;
				}
				fstring[i]=' ';
				float cx,cy,cz,xs,ys,zs;
				sscanf(fstring," 	%f#%f#%f#%f#%f#%f#",&cx,&cy,&cz,&xs,&ys,&zs); //get center_x center_y center_z size_x size_y size_z
				//printf(" 	%f	%f	%f	%f	%f	%f\n",cx,cy,cz,xs,ys,zs);
				
				char strn[50];
				strcat(box_data," --center_x ");
				snprintf(strn, sizeof(strn), "%.2f",cx);
				strcat(box_data,strn);
				strcat(box_data," --center_y ");
				snprintf(strn, sizeof(strn), "%.2f",cy);
				strcat(box_data,strn);
				strcat(box_data," --center_z ");
				snprintf(strn, sizeof(strn), "%.2f",cz);
				strcat(box_data,strn);	
				
				strcat(box_data," --size_x ");
				snprintf(strn, sizeof(strn), "%.2f",xs);				
				strcat(box_data,strn);			
				strcat(box_data," --size_y ");
				snprintf(strn, sizeof(strn), "%.2f",ys);				
				strcat(box_data,strn);				
				strcat(box_data," --size_z ");
				snprintf(strn, sizeof(strn), "%.2f",zs);				
				strcat(box_data,strn);
				break;
			}
			
		
		}
		
		
		
	}
	else
	{
		printf("\n Open file: ""./info_data/box.data"" failed \n");
		fclose(boxfp);
		exit(-1);
	}
	fclose(boxfp);
	return 0;
}









int covert2pdbqt(char *name)
{
	char command[300],cpname[100];
	command[0]='\0';
	strcat(command,"obabel ");
	strcat(command,name);
	int i=0;
	while(name[i]!='\0') 
	{
		cpname[i]=name[i];
		i++;
	}
	cpname[i]='\0';
	i=0;
	while(cpname[i]!='\0') i++; 
	while(cpname[i]!='.') i--;
	cpname[i]='\0';
	strcat(cpname,".pdbqt");
	strcat(command," -o pdbqt -O ");
	strcat(command,cpname);
	strcat(command," -b --partialcharge mmff94");
	system(command);
	return 0;
}



struct score_stats rundockf(char *receptor, char *ligand, int id) {
    char cpname[100];
    char command[300];
    struct score_stats stats = {0};
    int  DBL_MAX = 1000;
    // initialization stats
    stats.best_score = DBL_MAX;
    stats.worst_score = -DBL_MAX;
    stats.average_score = 0.0;
    stats.count = 0;
    
    // Construct ligand PDBQT file name
    int i = 0;
    while(ligand[i] != '\0') {
        cpname[i] = ligand[i];
        i++;
    }
    cpname[i] = '\0';
    i = 0;
    while(cpname[i] != '\0') i++; 
    while(cpname[i] != '.') i--;
    cpname[i] = '\0';	
    strcat(cpname, ".pdbqt");
    
    // Build Command
    command[0] = '\0';
    strcat(command, "vina --ligand ");
    strcat(command, cpname);
    
    char boxd[300];
    boxd[0] = '\0';
    get_boxf(receptor, boxd);	
    strcat(command, boxd);
    strcat(command, " --log  ");
    sprintf(cpname, "%d", id);
    strcat(command, cpname);
    strcat(command, ".log");
    
   // printf("Running docking: %s\n", command);
    
    
    // Hide Vina output
	#ifdef _WIN32
		strcat(command, " > nul 2>&1");  // Windows
	#else
    	strcat(command, " > /dev/null 2>&1");  // Linux/Mac
	#endif
        
    
    
    
    system(command);
    
    // Retrieve complete rating statistics from log files
    char log_filename[50];
    sprintf(log_filename, "%d.log", id);
    
    if (extract_scores(log_filename, &stats) != 0) {
        printf("Warning: Could not extract docking scores from %s\n", log_filename);
    }
    
    return stats;
}





void compute_descriptor(char *filename, struct atom_info *atom_info_data, struct hu_data *output_descriptor) {
    struct atoms *atom_data_local = NULL; 
    int atom_number_local, file_number_local;
    
    // Read file data
    get_data_from_file(filename, atom_info_data, &atom_data_local, &atom_number_local, &file_number_local);
    
    // Set descriptor name (extract from file name, remove path and extension)
    char *base_name = strrchr(filename, '/');
    if (base_name == NULL) base_name = filename;
    else base_name++;
    
    strcpy(output_descriptor->name, base_name);
    char *dot = strrchr(output_descriptor->name, '.');
    if (dot != NULL) *dot = '\0';
    
    //Calculate descriptors
    calculate(atom_data_local, atom_number_local, output_descriptor);
    
    // Release atomic data
    free(atom_data_local);
}














int split_xyz(char *file_name) {
    FILE *input_file = fopen(file_name, "r");
    if (input_file == NULL) {
        perror("Error opening input file");
        return -1;
    }

    char line[256];
    int file_count = 0;
    int atom_count = 0;
    int current_atoms = 0;
    FILE *output_file = NULL;

    while (fgets(line, sizeof(line), input_file)) {
        // Remove line breaks at the end of a line
        line[strcspn(line, "\n")] = 0;
        if (output_file == NULL) {
            atom_count = atoi(line);
            if (atom_count <= 0) {
                continue;
            }

            // Open a new output file
            file_count++;
            char output_name[50];
            snprintf(output_name, sizeof(output_name), "temp_%d.xyz", file_count);
            output_file = fopen(output_name, "w");
            if (output_file == NULL) {
                perror("Error creating output file\n");
                fclose(input_file);
                return -1;
            }

            // Write the number of atoms
            fprintf(output_file, "%d\n", atom_count);
            current_atoms = 0;
            if (fgets(line, sizeof(line), input_file)) {
                line[strcspn(line, "\n")] = 0;
                fprintf(output_file, "%s\n", line);
            } else {
                fclose(output_file);
                output_file = NULL;
                file_count--; 
                break;
            }
        } else {
            fprintf(output_file, "%s\n", line);
            current_atoms++;
            if (current_atoms >= atom_count) {
                fclose(output_file);
                output_file = NULL;
            }
        }
    }

    if (output_file != NULL) {
        fclose(output_file);
    }

    fclose(input_file);
    return file_count;
}




//  Preview Boltzmann distribution without filtering
// New function: Preview Boltzmann distribution without filtering
int preview_boltzmann_distribution(char *file_name) {
    int atom_number, file_number, count = 0, fcount = 0, linenumber = 0;
    FILE *fp;
    
    atom_number = get_xyz_atom_number(file_name);
    file_number = get_xyz_file_number(file_name);    
    
    if((atom_number <= 1) || (file_number <= 1)){
        printf("\n\nRead atom number %d or file number %d failed, file: %s\n", atom_number, file_number, file_name);
        return -1;
    }
    
    struct energy {
        int id;
        float energy;
    };
    
    struct energy *energy_idx;
    energy_idx = (struct energy *)malloc(file_number * sizeof(struct energy));
    char fs[300];
    
    fp = fopen(file_name, "r");
    if(fp == NULL){
        printf("\n\nCannot open file: %s\n", file_name);
        return -1;
    } else {
        while(fgets(fs, 255, fp) != NULL){
            if (strlen(fs) > 1 && (fs[strlen(fs)-2] == '\r')){
                fs[strlen(fs)-2] = '\n';
                fs[strlen(fs)-1] = '\0';
            }
            
            linenumber++;
            if((linenumber - ((atom_number + 2) * fcount)) == 2){
                float energy;
                int i = 0;
                while(fs[i] != ' ') {
                    fs[i] = ' ';
                    i++;
                }
                sscanf(fs, "%f", &energy);
                energy_idx[fcount].energy = energy;
                energy_idx[fcount].id = fcount + 1;
            }
            if((linenumber % (atom_number + 2)) == 0){
                fcount++;
            }
        }    
        fclose(fp);
    }
    
    // Sort by energy
    int tempid;
    float tempe;
    int i, j;
    for(i = 0; i <= file_number - 1; i++){
        for(j = i + 1; j <= file_number - 1; j++){
            if(energy_idx[i].energy >= energy_idx[j].energy){
                tempid = energy_idx[i].id;
                tempe = energy_idx[i].energy;
                energy_idx[i].id = energy_idx[j].id;
                energy_idx[i].energy = energy_idx[j].energy;
                energy_idx[j].id = tempid;
                energy_idx[j].energy = tempe;
            }
        }
    }
    
    float tempature = 298.15, R = 0.0019872, dis = 0;
    
    for(i = 0; i <= file_number - 1; i++) {
        dis = exp(-(energy_idx[i].energy - energy_idx[0].energy) / (tempature * R)) + dis;
    }
    
    // Display Boltzmann distribution preview
    
	printf("\n"
		" _______________________________________________________________________ \n"
		"|                                                                       |\n"
		"|                 === Boltzmann Distribution Preview ===                |\n"
		"|_______________________________________________________________________|\n\n");	    
    
	printf("      Conf    Energy (kcal/mol)    Dist (%%)    Cumulative (%%)\n");
	printf("_______________________________________________________________________\n");
    
    float cumulative_percent = 0.0;
    for(i = 0; i <= file_number - 1 && i < 20; i++) {  // Only show first 20
        float percent = exp(-(energy_idx[i].energy - energy_idx[0].energy) / (tempature * R)) / dis * 100;
        cumulative_percent += percent;
    printf("      %-4d       %12.3f         %8.3f        %8.3f\n", 
           energy_idx[i].id, 
           energy_idx[i].energy, 
           percent, 
           cumulative_percent);
    }
	printf("_______________________________________________________________________\n");    
    if(file_number > 20) {
        printf("... (Total %d conformations, showing first 20)\n", file_number);
        
        // Calculate total coverage
        float total_covered = 0.0;
        for(i = 0; i <= file_number - 1; i++) {
            total_covered += exp(-(energy_idx[i].energy - energy_idx[0].energy) / (tempature * R)) / dis * 100;
        }
        printf("Total distribution of all conformations: %.3f%%\n", total_covered);
    }
    
    // Provide recommended thresholds
    printf("\n=== Recommended Thresholds ===\n");
    printf("90%% coverage: ");
    cumulative_percent = 0.0;
    for(i = 0; i <= file_number - 1; i++) {
        float percent = exp(-(energy_idx[i].energy - energy_idx[0].energy) / (tempature * R)) / dis * 100;
        cumulative_percent += percent;
        if(cumulative_percent >= 90.0) {
            printf("Approx %.3f%% (select first %d conformations)\n", percent, i + 1);
            break;
        }
    }
    
    printf("95%% coverage: ");
    cumulative_percent = 0.0;
    for(i = 0; i <= file_number - 1; i++) {
        float percent = exp(-(energy_idx[i].energy - energy_idx[0].energy) / (tempature * R)) / dis * 100;
        cumulative_percent += percent;
        if(cumulative_percent >= 95.0) {
            printf("Approx %.3f%% (select first %d conformations)\n", percent, i + 1);
            break;
        }
    }
    
    printf("99%% coverage: ");
    cumulative_percent = 0.0;
    for(i = 0; i <= file_number - 1; i++) {
        float percent = exp(-(energy_idx[i].energy - energy_idx[0].energy) / (tempature * R)) / dis * 100;
        cumulative_percent += percent;
        if(cumulative_percent >= 99.0) {
            printf("Approx %.3f%% (select first %d conformations)\n", percent, i + 1);
            break;
        }
    }
    
    free(energy_idx);
    return 0;
}



float get_dis_with_preview(char *file_name) {
    char buffer_command[300];
    float bloz;
    
    // First show Boltzmann distribution preview
    if(preview_boltzmann_distribution(file_name) != 0) {
        printf("Cannot calculate Boltzmann distribution, using default value 0.5%%\n");
        return 0.5f;
    }
    
    printf("\n=== Boltzmann Distribution Filtering ===\n");
    printf("Enter the minimum distribution percentage threshold (0-100):\n");
    printf("Recommended: 0.5-5.0 (covers most important conformations)\n");
    printf("Input value: ");
    
    fgets(buffer_command, sizeof(buffer_command), stdin);
    sscanf(buffer_command, "%f", &bloz);
    
    // Validate input reasonableness
    if(bloz < 0.01f) {
        printf("Warning: Threshold %.3f%% is too low, may select too many conformations\n", bloz);
        printf("Confirm using this threshold? (y/n): ");
        fgets(buffer_command, sizeof(buffer_command), stdin);
        if(buffer_command[0] != 'y' && buffer_command[0] != 'Y') {
            return get_dis_with_preview(file_name); // Re-enter
        }
    } else if(bloz > 50.0f) {
        printf("Warning: Threshold %.3f%% is too high, may miss important conformations\n", bloz);
        printf("Confirm using this threshold? (y/n): ");
        fgets(buffer_command, sizeof(buffer_command), stdin);
        if(buffer_command[0] != 'y' && buffer_command[0] != 'Y') {
            return get_dis_with_preview(file_name); // Re-enter
        }
    }
    
    printf("Your selected Boltzmann distribution threshold: %.3f%%\n", bloz);
    return bloz;
}










int gen_conf(char *input_file) 
{
		int genc_method;
		char buffer_command[300],newname[100],newname1[100];
		float sel_bloz;
		int rotatornumb=0;
		rotatornumb=get_rotatable_bonds(input_file);
		
		
		
		
		printf("\n\n"
			" _______________________________________________________________________ \n"
			"|                                                                       |\n"
			"|                     === Conformation Generation ===                   |\n"
			"|_______________________________________________________________________|\n\n");	
	
		printf("1. Monte Carlo method (fast, may generate more invalid structures)\n");
		printf("2. Confab method (slower, generates fewer wrong structures)\n");
		printf("3. Skip conformation generation (use input structure directly)\n");
		printf("Please select method [1-3]: ");
	//	printf("There are two methods for generation conformation 1 Monte Carlo method 2 confab 3 None generation conformation\n");
		if(rotatornumb>=5) printf("The number of rotatable bonds in the molecule you submitted is %d. We strongly recommend that you do not perform conformational sampling and proceed with subsequent calculations directly (option [3]).\n",rotatornumb);
		if((0<=rotatornumb)&&(rotatornumb<5)) printf("Recommend using [2]confab, which runs slowly but generates fewer wrong stractures.\n");
		if(rotatornumb<0){
			printf("Get the number of rotatable bonds of %s failed,terminate calculation.",input_file);
			exit(-1);
		} 
		printf("Please enter the selected method:\n");		
		if(fgets(buffer_command, sizeof(buffer_command), stdin)==NULL){
			fprintf(stderr,"Error reading input.\n");
			return 1;
		}
		sscanf(buffer_command, "%d", &genc_method);	
		if(genc_method==1) printf("You chooed method is：[1] Monte Carlo method \n");	
		else if(genc_method==2) printf("You chooed method is：[2] confab method \n");
		else if(genc_method==3) printf("You chooed method is：[3] None generation conformation \n");
		else{
		//	printf("You chooed method is：[%], undefined exit\n");
			fprintf(stderr, "You chooed method is：[%d], undefined exit\n",genc_method);
			return -1;
		}
		if((genc_method==1)||(genc_method==2))
		{
			
		//	printf(" MMFF94 force field for conformational optimization\n");

			if(genc_method==1){
	 			cal_conformer(input_file);			
			}
			if(genc_method==2){	
				if(cal_confab(input_file)==1) return 1;	
			}
			add_prefix(newname, input_file, "cos_");
			chang_suffix2xyz(newname);
		
			cal_obabel_mmf(newname);
 		//	sel_bloz= get_dis();
 			
			add_prefix(newname1, newname, "mm_");
			
			sel_bloz = get_dis_with_preview(newname1);
			if(get_energy(newname1,sel_bloz)==1) return 1;	
			
			int filenumber=0;
			int i=0;
			while(newname1[i]!='\0')i++;
			while(newname1[i]!='.')i--;
			newname1[i]='\0';
			strcat(newname1,"_choosed.xyz");
			char command[300];
			command[0]='\0';
			strcat(command,"obabel -i xyz  ");
			strcat(command,newname1);

			i=0;
			while(newname1[i]!='\0')i++;
			while(newname1[i]!='.')i--;
			newname1[i]='\0';
			strcat(newname1,"_dh.xyz");
			strcat(command," -o xyz -O ");
			strcat(command,newname1);
			strcat(command," -d");
			system(command);
			filenumber=split_xyz(newname1);
			return filenumber;			
			
		}
		if(genc_method==3) return 1;
	
}
	





int cal_dis(struct atoms *atom_data,int atom_number,float x,float y,float z)
{

int i,fit=0 ,index_min;
float distance,min_dis;

for(i=0;i<=atom_number-1;i++)
{
	if(atom_data[i].ele>1)
	{
		min_dis=sqrt(  (atom_data[i].xa-x)*(atom_data[i].xa-x) + (atom_data[i].ya-y)*(atom_data[i].ya-y)  +  (atom_data[i].za- z)*(atom_data[i].za-z));
		break;
	}
}

//min_dis=sqrt(  (atom_data[0].xa-x)*(atom_data[0].xa-x) + (atom_data[0].ya-y)*(atom_data[0].ya-y)  +  (atom_data[0].za- z)*(atom_data[0].za-z));
index_min=0;
for(i=0;i<=atom_number-1;i++)
{	
	if(atom_data[i].ele>1)
	{
		distance=sqrt(  (atom_data[i].xa-x)*(atom_data[i].xa-x) + (atom_data[i].ya-y)*(atom_data[i].ya-y)  +  (atom_data[i].za- z)*(atom_data[i].za-z));
		if(min_dis>distance){
			min_dis=distance;
			index_min=i;	
		} 
		if(distance<=atom_data[i].rvdw) fit++;	
	}

}

if(fit==0){
	if(min_dis<=atom_data[index_min].rvdw+2.5)
	return 0;
	else return 100;
} 
if(fit==1) return 1;
if(fit>1) return  2;
else return -1;
}









struct min_max  get_min_max(struct atoms *atom_data,int atom_number)
{
	struct min_max minmax;
	int i;
	float xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
	xmin=atom_data[0].xa;xmax=atom_data[0].xa;
	ymin=atom_data[0].ya;ymax=atom_data[0].ya;
	zmin=atom_data[0].za;zmax=atom_data[0].za;
	for(i=0;i<=atom_number-1;i++)
	{
		xmin= (atom_data[i].xa < xmin) ? atom_data[i].xa : xmin;
		xmax= (atom_data[i].xa > xmax) ? atom_data[i].xa : xmax;
		ymin= (atom_data[i].ya < ymin) ? atom_data[i].ya : ymin;
		ymax= (atom_data[i].ya > ymax) ? atom_data[i].ya : ymax;
		zmin= (atom_data[i].za < zmin) ? atom_data[i].za : zmin;
		zmax= (atom_data[i].za > zmax) ? atom_data[i].za : zmax;		
	}
	minmax.xmin=xmin;
	minmax.xmax=xmax;
	minmax.ymin=ymin;
	minmax.ymax=ymax;
	minmax.zmin=zmin;
	minmax.zmax=zmax;
	return minmax;
}




int cal_adjacent(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number)
{
	int xadd1,xsub1,yadd1,ysub1,zadd1,zsub1,adja=0,count=0,gen_number1=0,gen_number2=0;
	if(grid_data[number].state==0)
	{
		
		//  x-x
		xadd1=number+(y_number*z_number);
		xsub1=number-(y_number*z_number);
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
		


		//y-y
		yadd1=number+z_number;
		ysub1=number-z_number;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0)&&(grid_data[yadd1].state!=100)&&(grid_data[ysub1].state!=100))
			count++;
		}
		
	

		
		//z-z
		zadd1=number+1;
		zsub1=number-1;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0)&&(grid_data[zadd1].state!=100)&&(grid_data[zsub1].state!=100))
			count++;
		}	
		
		//x-1,y+1,z---------x+1,y-1,z
		xadd1=number+(y_number*z_number)-z_number;
		xsub1=number-(y_number*z_number)+z_number;		
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
		
		
				
		// x-1,y-1,z---------x+1,y+1,z
		xadd1=number+(y_number*z_number)+z_number;
		xsub1=number-(y_number*z_number)-z_number;		
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}	
		
		
			
		// x+1,y+1,z+1--------	x-1,y-1,z-1	
		xadd1=number+(y_number*z_number)+z_number+1;
		xsub1=number-(y_number*z_number)-z_number+1;
		
	
					
		if(    ((xadd1>0)&&(xadd1<=grid_number-1)) &&
		       ((xsub1>0)&&(xsub1<=grid_number-1))
		  )
		{
		
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				
				count++;
			}
			
		}		
	

		// x,y+1,z+1--------	x,y-1,z-1	
		xadd1=number+z_number+1;
		xsub1=number-z_number-1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}		

		// x-1,y+1,z+1--------	x+1,y-1,z-1	
		xadd1=number+(y_number*z_number)-z_number-1;
		xsub1=number-(y_number*z_number)+z_number+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}		

		// x-1,y,z+1--------	x+1,y,z-1	
		xadd1=number+(y_number*z_number)-1;
		xsub1=number-(y_number*z_number)+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}	

		// x-1,y-1,z+1--------	x+1,y+1,z-1	
		xadd1=number+(y_number*z_number)+z_number-1;
		xsub1=number-(y_number*z_number)-z_number+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}

		// x,y-1,z+1--------	x,y+1,z-1	
		xadd1=number+z_number-1;
		xsub1=number-z_number+1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}


		// x+1,y-1,z+1--------	x-1,y+1,z-1	
		xadd1=number+(y_number*z_number)-z_number+1;
		xsub1=number-(y_number*z_number)+z_number-1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}

		// x+1,y,z+1--------	x-1,y,z-1	
		xadd1=number+(y_number*z_number)+1;
		xsub1=number-(y_number*z_number)-1;			
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}


		if(count>=1) return 3;
		else return 1;	
	}
	else return 0;

}


int cal_adjacent1(int number, struct grid *grid_data,int grid_number,int x_number,int y_number,int z_number)
{
	int xadd1,xsub1,yadd1,ysub1,zadd1,zsub1,adja=0,count=0,gen_number1=0,gen_number2=0;
	if(grid_data[number].state==0)
	{
		//  x-x
		xadd1=number+(y_number*z_number);
		xsub1=number-(y_number*z_number);
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
		
		xadd1=number+(y_number*z_number)*1;
		xsub1=number-(y_number*z_number)*2;
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}
			
		
		xadd1=number+(y_number*z_number)*2;
		xsub1=number-(y_number*z_number)*1;
		if(((xadd1>0)&&(xadd1<=grid_number-1)) &&((xsub1>0)&&(xsub1<=grid_number-1)))
		{
			if((grid_data[xadd1].state!=0)&&(grid_data[xsub1].state!=0) )
			{
				count++;
			}
			
		}

		//y-y
		yadd1=number+z_number;
		ysub1=number-z_number;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}
		
		yadd1=number+z_number*2;
		ysub1=number-z_number*2;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}		

		yadd1=number+z_number*1;
		ysub1=number-z_number*2;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}		
				
		yadd1=number+z_number*2;
		ysub1=number-z_number*1;
		if(((yadd1>0)&&(yadd1<=grid_number-1)) &&((ysub1>0)&&(ysub1<=grid_number-1)))
		{
			if((grid_data[yadd1].state!=0)&&(grid_data[ysub1].state!=0) )
			count++;
		}				
		//z-z
		zadd1=number+1;
		zsub1=number-1;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}	
		zadd1=number+2;
		zsub1=number-2;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}	
		zadd1=number+1;
		zsub1=number-2;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}	
		
						
		zadd1=number+2;
		zsub1=number-1;		
		
		if(((zadd1>0)&&(zadd1<=grid_number-1)) &&((zsub1>0)&&(zsub1<=grid_number-1)))
		{
			if((grid_data[zadd1].state!=0)&&(grid_data[zsub1].state!=0) )
			count++;
		}
		if(count>=1)
			return 3;
		else return 0;	
	}
	else return 0;
}








int file_type(char *file_name)
{
	char *mol2=".MOL2";
	char *sdf=".SDF";
	char *pdb=".PDB";
	char *xyz=".XYZ";
	char *p_name;
	char suffix_name[6];
	int i=0;
	p_name=file_name;
	while(*p_name!='\0') p_name++;
	while(*p_name!='.') p_name--;
	while(*p_name!='\0')
	{	
		if(i>=5)
		{	
			printf("\nFile format error\n");
			return -1;
		 } 
		suffix_name[i]=toupper(*p_name);
		p_name++;
		i++;
	} 
	suffix_name[i]='\0';
	if(strcmp(suffix_name, mol2)==0) return 1; 
	else if(strcmp(suffix_name, pdb)==0) return 2; 
	else if(strcmp(suffix_name, sdf)==0) return 3; 
	else if(strcmp(suffix_name, xyz)==0) return 4;
	else
	{
			printf("\nFile format error\n");
			return -1;		
	 } 

}


void usage()
{   
	printf("calmvs [-i<in file name>] [-b < center x, center y, center z, size in the x dimension (Angstroms), size in the y dimension (Angstroms) size in the z dimension (Angstroms)>][Options] [-t<calculate Polar Surface Area>][Options] [-o<out file name>][Options]");
	printf("Try  -H option for more information.");
}
void help()
{   
	printf("calmvs calculates molecule surface volume and compares molecular similarity\n\n\n");
	printf("Usage\n");
	printf("calmvs [-i<in file name>] [-b < center x, center y, center z, size in the x dimension (Angstroms), size in the y dimension (Angstroms) size in the z dimension (Angstroms)>][Options] [-t<calculate Polar Surface Area>][Options] [-o<out file name>][Options]");
	printf("The extension of a file decides the format,Only supports mol2, pdb, sdf formats\n");	
}


int get_atom_info(struct atom_info *atom_info_data)
{
	FILE *fp;
	char fstring[256];
	int atom_number,line_number=0,i,j;
	char *name="./info_data/ele_data.txt";
	fp=fopen(name,"r");
	
	if(fp!=NULL){
		while(fgets(fstring,255,fp)!=NULL){
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    		}
			if(fstring[0]=='#') ;
			else{
				int index; 
				char atom[5],rvdw[10],mw[10];
				sscanf(fstring,"%d %s %s %s",&index,atom,rvdw,mw);
    			char *endptr;
    			float frvdw,fmw;
     			frvdw= strtof(rvdw, &endptr);
				if (*endptr != '\0') { // Check if the entire string has been converted
        			printf("Conversion failed, unconverted parts: %s\n", endptr);
    			} else {
        		;	
    			}
     			fmw= strtof(mw, &endptr);
				if (*endptr != '\0') { // Check if the entire string has been converted
        			printf("Conversion failed, unconverted parts: %s\n", endptr);
    			} else {
        			;
    			}
    			if((0<index)&&(index<199)){
    				atom_info_data[index].falg=1;
    				atom_info_data[index].id=index;
    				atom_info_data[index].mw=fmw;
    				atom_info_data[index].vdwr=frvdw;
    				int ia=0;
    				for(ia=0;ia<=4;ia++) atom_info_data[index].atoms[ia]=atom[ia];
    		
				}else{
					printf("Get atom number data from ele_data.txt failed,  atom number<0 or atom number > 300 \n");
				}
			
			}
		}
					
	}
	else{
		printf ("/n Open file %s failed, please check it\n",name);
		exit(-1);
	}
	return 0;
}


int get_data(struct atoms **atom_data,struct atom_info *atom_info_data,int count, float x,float y,float z,char *ele_s)
{
	
	(*atom_data)[count].xa=x;
	(*atom_data)[count].ya=y;
	(*atom_data)[count].za=z;
	
	if((ele_s[1])=='\0'){
		int temp_i=0,flag=0;
		for(temp_i=0;temp_i<=199;temp_i++){
			if(atom_info_data[temp_i].falg==1){
				if((ele_s[0]==atom_info_data[temp_i].atoms[0])&&(atom_info_data[temp_i].atoms[1]=='\0')){
					(*atom_data)[count].ele=atom_info_data[temp_i].id;
					(*atom_data)[count].rvdw = atom_info_data[temp_i].vdwr;
					flag++;
					break;}else;											
				}else ; 
		}
		if(flag==0){
			printf("\n Error cant get atom type #%s#\n",ele_s);
			return -1;}			
	}else if(ele_s[2]=='\0'){
		int temp_i=0,flag=0 ;
		for(temp_i=0;temp_i<=199;temp_i++){
			if(atom_info_data[temp_i].falg==1){
				if((ele_s[0]==atom_info_data[temp_i].atoms[0])&&(atom_info_data[temp_i].atoms[1]==ele_s[1])){
					(*atom_data)[count].ele=atom_info_data[temp_i].id;
					(*atom_data)[count].rvdw = atom_info_data[temp_i].vdwr;
					flag++;
					break;}												
			}
		}
		if(flag==0){
			printf("\n Error cant get atom type %s\n",ele_s);
			return -1;}									
	}else{
			printf("\nError wrong atom type\n");
			return -1;									
	}	
	
	return 0;
}











int get_data_from_file(char *file_name,struct atom_info *atom_info_data,struct atoms **atom_data,int *p_atom_number, int *xyz_file_number)
{
	FILE *fp;
	int atom_number=0,i=0,j=0,count=0,file_number=0;
	float x,y,z;
	char  fstring[356], ele_s[8],cotain[300],ele[8], ele2[8],atom_arry[5];
	*xyz_file_number=1;
	fp=fopen( file_name ,"r");
	
	if(fp==NULL){
			printf ("/n Open file %s failed\n", file_name );
			return -1;		
	}
	
	if(file_type( file_name )==-1){
		printf ("/n Not supporting file format of file %s \n", file_name );
		return -1;
	}	
	if(file_type( file_name )==1){
		int start=0,end=0,lock_mol2=0,number_mol2,number_mol2_1;
		float charge;
		char new_line[80],mol2_name[15];
		char *atom_start="@<TRIPOS>ATOM",*atom_end="@<TRIPOS>BOND",*linep,*linep1;
		
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if(strncmp(fstring,atom_start,13)==0) start=1;
				if(strncmp(fstring,atom_end,13)==0) end=1;
				if((start==1)&&(end==0)&&(strncmp(fstring,"@",1)!=0)){	
					atom_number++;
				}
			}
		}else{
			printf ("/n Open file %s failed ", file_name );
			fclose(fp);
			return -1;		
		}
		rewind(fp);
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
		else {
			printf ("/n ERROR read the count of atoms failed \n");
			return -1;
		}
		start=0;
		end  =0;
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if(strncmp(fstring,atom_start,13)==0) start=1;
				if(strncmp(fstring,atom_end,13)==0) end=1;
				if((start==1)&&(end==0)&&(strncmp(fstring,"@",1)!=0)){	
              	  linep=fstring;
             	   i=0;
            		while(*linep!='\0'){                  	
                		if(isgraph(*linep)!=0){
							   if(lock_mol2==0){
					   				linep1=linep;
                   					while(isgraph(*linep1)!=0){
                    					cotain[i]=*linep1;
                    					linep1++;
                   						i++;                  				
									}
									cotain[i]=' ';
									i++;						   	
							   }	
                   				lock_mol2=1;
						}
						if(isgraph(*linep)==0) lock_mol2=0;
						linep++;
					}
					cotain[i]='\0';
              		sscanf(cotain,"%d %s %f %f %f %s %d %s %f",&number_mol2,ele,&x,&y,&z,ele2,&number_mol2_1,mol2_name,&charge);
               		i=0;
              		while(ele2[i]!='\0'){
                  	 		ele2[i]=toupper(ele2[i]);
                  	 		if(	ele2[i]=='.'){
                   			ele2[i]='\0';
                   			break;
					  	 } else;
                   		i++;
					}
              		if(count<atom_number){
						get_data(atom_data,atom_info_data,count,x,y,z,ele2);	    
						count++;							
					}else{
						printf("\nError read atom data error\n");
						exit(-1);
					}		
				}
			}
		}else{
			printf ("/n Open file %s failed ", file_name );
			fclose(fp);
			return -1;		
		}
		
		*p_atom_number=atom_number;
	} else if(file_type( file_name )==2)
	{
		char *atom="ATOM",*hetatom="HETATM";	
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if((strncmp(fstring,atom,4)==0)||(strncmp(fstring,hetatom,6)==0)){
					atom_number++;
				}
			}		
		}else{
		printf ("/n Open file %s failed ", file_name );
		fclose(fp);
		return -1;
		}
		rewind(fp);
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
		else {
			printf ("/n ERROR read the count of atoms failed \n");
			return -1;
		}
		if(fp!=NULL){
			while ((fgets(fstring,119,fp))!=NULL){
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if((strncmp(fstring,atom,4)==0)||(strncmp(fstring,hetatom,6)==0)){
					if(fstring[76]==' '){
						ele_s[0]=toupper(fstring[77]);
						ele_s[1]='\0';
					//	printf("%s\n",ele_s);
					}else if(((fstring[76]>='A')&&(fstring[76]<='Z'))||((fstring[76]>='a')&&(fstring[76]<='z'))){
						ele_s[0]=toupper(fstring[76]);
						if(fstring[77]==' ')ele_s[1]='\0';
						else if(((fstring[77]>='A')&&(fstring[77]<='Z'))||((fstring[77]>='a')&&(fstring[77]<='z'))){
							ele_s[1]=toupper(fstring[77]);
							ele_s[2]='\0';							
						}else ele_s[1]='\0';

					//	printf("%s\n",ele_s);						
					}else{
							printf("\nError wrong PDB file data colume 77 and colume 78 no elements data\n");
							return -1;						
					}	
					for(i=0;i<=29;i++) fstring[i]=' ';
					for(i=54;i<=79;i++) fstring[i]=' ';
				//	printf("#%s#\n",ele_s);
					sscanf(fstring,"                                %f  %f  %f                        ",&x,&y,&z);
					if(count<atom_number){
							get_data(atom_data,atom_info_data,count,x,y,z,ele_s);
							count++;
					}
					else{
						printf("\nError read atom data error\n");
						exit(-1);
					}
				}
			}		
		}
		else{
			printf ("/n Open file %s failed \n", file_name );
			fclose(fp);
			return -1;
		}
		
		*p_atom_number=atom_number;
	}else if (file_type( file_name )==3)
	{
		int line_count=0;	
		count=0; 
		if(fp==NULL){
			printf("Open file failed\n");
			return -1;
		}
		while(fgets(fstring,255,fp)!=NULL){  
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		} 
			line_count++;	
			if(line_count==4){	
			fstring[3]='\0';
			sscanf(fstring,"%d",&atom_number);
			break;
			} 
		}
		*p_atom_number=atom_number;
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*sizeof(struct atoms ));
		else {
			printf ("\n ERROR read the count of atoms failed atom number %d\n",atom_number);
			return -1;
		}
		line_count=1;
		while(fgets(fstring,255,fp)!=NULL){  
		 	if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		}
			if(line_count<=atom_number)	{	
				i=0;
				while(isalpha(fstring[i])==0) i++;
				while(isalpha(fstring[i])) i++;
				fstring[i]='\0';
				i=0; 
				j=0;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]=' ';
				j++;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]=' ';
				j++;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]=' ';
				j++;
				while(isspace(fstring[i])) i++;
				while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
					cotain[j]=fstring[i];
					i++;
					j++;
				}
				cotain[j]='\0';	
				sscanf(cotain,"%f %f %f %s",&x,&y,&z,ele);
				i=0;
           	 while(ele[i]!='\0'){
            	    ele[i]=toupper(ele[i]);
            	    i++;
				}
				if(count<atom_number){
					get_data(atom_data,atom_info_data,count,x,y,z,ele);	    
					count++;							
				}else {
					printf("\nError read atom data error\n");
					exit(-1);
				}		
			}	else break;
			line_count++;		
		}
		
		*p_atom_number=atom_number;
	}	else if (file_type( file_name )==4)
	{
		
		int line_count=0; 
		while(fgets(fstring,255,fp)!=NULL){
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		}  	
			line_count++;
			if(line_count==1){
				i=0;
				j=0;
				while(isspace(fstring[i])) i++;
				while(isspace(fstring[i])==0){
					atom_arry[j] =fstring[i];
					i++;
					j++;					
				}
				atom_arry[j]='\0';
				sscanf(atom_arry,"%d",&atom_number);
				
			}
			if(line_count==2) break;
		}
		rewind(fp);
		line_count=0;
		while(fgets(fstring,255,fp)!=NULL){
			if(strlen(fstring)>=2) line_count++;		
		}
		file_number=line_count/(atom_number+2);
		*p_atom_number=atom_number;
		*xyz_file_number=file_number;	
		count=0; 
		rewind(fp);
		if((atom_number>0)&&(atom_number<=100000)) *atom_data=(struct atoms *)malloc(atom_number*file_number*sizeof(struct atoms ));
		else {
			printf ("\n ERROR read the count of atoms failed \n");
			return -1;
		}	
		line_count=0;
		line_count=atom_number+2;
		count=0;
		
		
		while(fgets(fstring,255,fp)!=NULL){
			if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      			fstring[strlen(fstring)-2]='\n';
      			fstring[strlen(fstring)-1]='\0';
    		}
		   
			line_count++;
			if(((line_count-1)%(atom_number+2)!=0)&&((line_count-2)%(atom_number+2)!=0)){
				if(line_count<=(atom_number+2)*(file_number+1)){
					i=0; 
					j=0;				
					if(fstring[0]==' ') while(isspace(fstring[i])) i++;
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						if(isalpha(fstring[i])){
							cotain[j]=fstring[i];
					//	printf("cotain[%d] --[%c]\n",j,cotain[j]);
							i++;       
							j++;							
						}
						else i++;
					}
					cotain[j]=' ';
					j++;
					if(fstring[i]=' ')	while(isspace(fstring[i])) i++;
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						cotain[j]=fstring[i];
						i++;  //x
						j++;
					}
					cotain[j]=' ';
					j++;
					
					while(isspace(fstring[i])) i++;
					
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						cotain[j]=fstring[i];
						i++;   //y
						j++;
					}
					cotain[j]=' ';
					j++;
					while(isspace(fstring[i])) i++;
					while((isspace(fstring[i])==0)&&(fstring[i]!='\0')){
						cotain[j]=fstring[i];
						i++;
						j++;
					}
					j++;
					cotain[j]='\0';	
					sscanf(cotain,"%s %f %f %f",ele,&x,&y,&z);
					i=0;
           			while(ele[i]!='\0'){
              			ele[i]=toupper(ele[i]);
               			i++;
					}
					ele[2]='\0';
					if(count<=atom_number*file_number-1){
						get_data(atom_data,atom_info_data,count,x,y,z,ele);	    
						count++;					
					}	
				}else break;			
			}	
		}
		
		
	}else{
		printf("\nError file type error\n");
		fclose(fp);
		exit(-1);		
	}
	fclose(fp);	
	return 	0;		
}


int calculate(struct atoms *atom_data, int atom_number,struct hu_data *hudata)
{
	float grid_gap;	
	long long int i,j ,k,l, box_x,box_y,box_z,count=0,s_cal_dis2,cal_n,cal_nadd1,cal_nsub1;
	float x,y,z,box_vol;
	int grid_n,line_count=0,*noblank_index;
	struct data_index *data_index;
	int box_a,box_b,box_c,temp;
	struct min_max minmax;
	struct timespec start, end;


	minmax= get_min_max(atom_data,atom_number);
    grid_gap=0.3800;

	box_x=(int)((minmax.xmax-minmax.xmin+16.0)/grid_gap);
	box_y=(int)((minmax.ymax-minmax.ymin+16.0)/grid_gap);
	box_z=(int)((minmax.zmax-minmax.zmin+16.0)/grid_gap);
	grid_data=(struct grid *)malloc(box_x*box_y*box_z*sizeof(struct grid));
	if ( grid_data == NULL) {
		fprintf(stderr, "memory allocation failed\n");
    	exit(EXIT_FAILURE);
	}	
  
	data_index=(struct data_index *)malloc(box_x*box_y*box_z*sizeof(struct data_index));
	if (data_index == NULL) {
		fprintf(stderr, "memory allocation failed\n");
    	exit(EXIT_FAILURE);
	}

	//printf("xmin%f xmax %f\n ymin %f ymax %f\n zmin %f zmax%f\n",minmax.xmin,minmax.xmax, minmax.ymin,minmax.ymax, minmax.zmin,minmax.zmax);

	count=0;
	x=minmax.xmin-8.0;
	for(i=0;i<=(box_x-1);i++){	
		y=minmax.ymin-8.0;
		for(j=0;j<=(box_y-1);j++){	
			z=minmax.zmin-8.0;
			for(k=0;k<=(box_z-1);k++){	
				if(count<=(box_x*box_y*box_z))	{
					grid_data[count].x=x;
					grid_data[count].y=y;
					grid_data[count].z=z;
				}
				else break;
				count++;
				z=z+grid_gap;
			}
			y=y+grid_gap;
		}	
		x=x+grid_gap;
	}	
	count=(box_x)*(box_y)*(box_z);
	for(i=0;i<=count-1;i++) data_index[i].wsp=0;


for(i=0;i<=count-1;i++)
{
grid_data[i].state   =  cal_dis(atom_data,atom_number,grid_data[i].x,grid_data[i].y,grid_data[i].z);		
if(grid_data[i].state==-1){
	fprintf(stderr,"\nError unsupport atom type.\n");
	return -1;

}	
else {
	if(grid_data[i].state<=2){
		data_index[i].wsp=1;
		data_index[i].x=grid_data[i].x;
		data_index[i].y=grid_data[i].y;
		data_index[i].z=grid_data[i].z;	
	} 
	}
}

int number_of_no_blank=0;
for(i=0;i<=count-1;i++){
	if(data_index[i].wsp==1) number_of_no_blank++;
}

noblank_index=(int *)malloc( number_of_no_blank*sizeof(int));
if (noblank_index == NULL) {
	fprintf(stderr, "memory allocation failed\n");
	exit(EXIT_FAILURE);
}

j=0;
for(i=0;i<=count-1;i++){
	if(data_index[i].wsp==1){
		noblank_index[j]=i;
		j++;
	}
}
for(i=0;i<=count-1;i++)
{

if(grid_data[i].state==100){
grid_data[i].state=0;

}

}



for(i=0;i<=number_of_no_blank-1;i++)
{

grid_data[noblank_index[i]].psa =0;
}




for(i=0;i<=number_of_no_blank-1;i++)
{
//	printf("%d - %d\n",i,cal_adjacent(i, grid_data,(box_x*box_y*box_z),box_x,box_y,box_z));
	if (cal_adjacent(noblank_index[i], grid_data,(box_x*box_y*box_z)-1,box_x,box_y,box_z)==3) 
	{    
		grid_data[noblank_index[i]].state=2;
	
	}
}












for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].state==1) 
	{    
		grid_data[noblank_index[i]].state=2;
	
	}
}
// x
for(i=0;i<=box_z-1;i++)
{
	for(j=0;j<=box_y-1;j++)
	{
		for(k=0;k<=box_x-1;k++)
		{
			cal_n=(k*box_y*box_z)+(j*box_z)+i;
			cal_nadd1=((k+1)*box_y*box_z)+(j*box_z)+i;
			cal_nsub1=((k-1)*box_y*box_z)+(j*box_z)+i;
			if((cal_nsub1>=0)&&(cal_nsub1<=(box_y*box_z*box_x-1)))
			{
				if((cal_nadd1>=0)&&(cal_nadd1<=(box_y*box_z*box_x-1)))
				{
					if((cal_n>=0)&&(cal_n<=(box_y*box_z*box_x-1)))
					{
						if(grid_data[cal_nadd1].state!=1)
						{
							if(grid_data[cal_nsub1].state!=1)
							{
							if(grid_data[cal_nadd1].state != grid_data[cal_nsub1].state ){
								grid_data[cal_n].state=1;
							} 
							
							}
						}
						

					}
				}
			}
			
			
			
			
			
		
		}
	  
	}  
		

}

for(i=0;i<=box_x-1;i++)
{
	for(j=0;j<=box_z-1;j++)
	{
		for(k=0;k<=box_y-1;k++)
		{
			cal_n=(i*box_y*box_z)+(k*box_z)+j;
			
	
			
			
			
			cal_nadd1=(i*box_y*box_z)+((k+1)*box_z)+j;
			cal_nsub1=(i*box_y*box_z)+((k-1)*box_z)+j;

				
			
			
			if((cal_nsub1>=0)&&(cal_nsub1<=(box_y*box_z*box_x-1)))
			{
				if((cal_nadd1>=0)&&(cal_nadd1<=(box_y*box_z*box_x-1)))
				{
					if((cal_n>=0)&&(cal_n<=(box_y*box_z*box_x-1)))
					{
						if(grid_data[cal_nadd1].state!=1)
						{
							if(grid_data[cal_nsub1].state!=1)
							{
								if(grid_data[cal_nadd1].state != grid_data[cal_nsub1].state ) 							
							 		 grid_data[cal_n].state=1;
							}
						}
						

					}
				}
			}
			
			
			
			
			
		
		}
	  
	}  
		

}

//z

for(i=0;i<=box_x-1;i++)
{
	for(j=0;j<=box_y-1;j++)
	{
		for(k=0;k<=box_z-1;k++)
		{
			cal_n=    (i*box_y*box_z)+(j*box_z)+k;
			cal_nadd1=(i*box_y*box_z)+(j*box_z)+k+1;
			cal_nsub1=(i*box_y*box_z)+(j*box_z)+k-1;
	
			if((cal_nsub1>=0)&&(cal_nsub1<=(box_y*box_z*box_x-1)))
			{
				if((cal_nadd1>=0)&&(cal_nadd1<=(box_y*box_z*box_x-1)))
				{
					if((cal_n>=0)&&(cal_n<=(box_y*box_z*box_x-1)))
					{
						if(grid_data[cal_nadd1].state!=1)
						{
							if(grid_data[cal_nsub1].state!=1)
							{
								if(grid_data[cal_nadd1].state != grid_data[cal_nsub1].state ) 
									grid_data[cal_n].state=1;
							}
						}
						

					}
				}
			}
			
			
			
			
			
		
		}
	  
	}  
		

}




count =0;
for(i=0;i<=number_of_no_blank-1;i++)
{
	
	if (grid_data[noblank_index[i]].state==1) 
	{    
		count++;
	
	}
}
float suface;
suface=count*grid_gap;
suface=suface*grid_gap;


suface = 1.26100 * suface + 13.50600;



hudata->v1=suface;




count=0;
for(i=0;i<=number_of_no_blank-1;i++)
{
if((grid_data[noblank_index[i]].state==2)|| (grid_data[noblank_index[i]].state==1))
	{    
		count++;	
	}
}

float volume;
volume=count*grid_gap;
volume=volume*grid_gap;
volume=volume*grid_gap;
hudata->v2=volume;
cal_moment( hudata ,grid_data, box_x*box_y*box_z , grid_gap, volume, suface);



free(grid_data);

 free(data_index );
 free(noblank_index);

return 0;	
}

int is_obabel_installed() {
    int ret;
#ifdef _WIN32
    // Windows uses the 'where' command to find obabel.exe
    ret = system("where obabel.exe > nul 2>&1");
#else
    // Unix/Linux/MacOS use which command to find obabel
    ret = system("which obabel > /dev/null 2>&1");
#endif
    return (ret == 0);
}

int get_xyz_atom_number(char *name)
{
	FILE *fp;
	char string[300];
	int file_numb;
	fp=fopen(name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",name);
		return -1;
	}
	else{
		if(fgets(string,255,fp)!=NULL){
			if (strlen(string) > 1 && (string[strlen(string)-2]=='\r')){
      			string[strlen(string)-2]='\n';
      			string[strlen(string)-1]='\0';
    		}
			sscanf(string,"%d",&file_numb);
			fclose(fp);
			return file_numb;
		}
		else{
			fclose(fp);
			return -1;
		} 
	}
}

int get_xyz_file_number(char *name)
{
	FILE *fp;
	char string[300];
	int atom_numb=0,linenumber=0;
	atom_numb=get_xyz_atom_number(name);
	if(atom_numb==-1){
		printf("\n\n get %s atom number failed\n",name);
		return -1;
	}
	fp=fopen(name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",name);
		return -1;
	}
	else{
		while(fgets(string,255,fp)!=NULL){
 			linenumber++;
		}
		fclose(fp);
	}
	return linenumber/(atom_numb+2);
}







int cal_confab(char *file_name)
{
	if(is_obabel_installed()!=1){
		printf("The presence of obabel was not detected on this device.\n Please install, recommended version: Open Babel 3.0.0\n");
		exit( -1);	
	}
	int conf=1000,i,conf_number=0,total_conf=0,round=0, rotbond=0;
	float rcutoff=0.01,ecutoff=5.0;
	char call_confab[200],string[200],*gen="..generated";
	char number[20],out_name[50];
	
	out_name[0]='\0';
	strcat(out_name,"cos_");
	strcat(out_name,file_name);
	i=0;
	while(out_name[i]!='\0') i++;
	while(out_name[i]!='.') i--;
	out_name[i]='\0';
	strcat(out_name,".xyz");	
	
	remove(out_name);
	remove("confab.log");	
	
	
	
	FILE *log_fp;
	for(rcutoff=0.5;rcutoff>=0.000;rcutoff=rcutoff-0.1) 
	{

			ecutoff=50.0;
		
			round++;
			call_confab[0]='\0';
			strcat(call_confab,"obabel ");
			strcat(call_confab,file_name);
			strcat(call_confab,"  -O  ");
			strcat(call_confab,out_name);
			strcat(call_confab," --confab --verbose ");   //--conf 1000000
			sprintf(number," %.3f ",rcutoff);
			strcat(call_confab," --rcutoff ");
			strcat(call_confab,number);	
			sprintf(number," %.3f ",ecutoff);		
			strcat(call_confab," --ecutoff ");
			strcat(call_confab,number);
			strcat(call_confab,">> confab.log");
			printf("%s\n",call_confab);
			system(call_confab);
			log_fp=fopen("confab.log","r");
			
			if(rcutoff<=0.1)
			{
				
				fclose(log_fp);
				return 0;				
			}
			if(log_fp!=NULL)
			{
				while(fgets(string,199,log_fp)!=NULL)
				{
					if (strlen(string) > 1 && (string[strlen(string)-2]=='\r')){
      					string[strlen(string)-2]='\n';
      					string[strlen(string)-1]='\0';
    				}
					
					
					if(strncmp(string,"..number of rotatable bonds =",29)==0)
					{
						sscanf(string,"..number of rotatable bonds = %d",&rotbond);
					}
					
					
					if(strncmp(string,"..tot conformations =",21)==0)
					{
						sscanf(string,"..tot conformations = %d",&total_conf);	
					}
					
					if(strncmp(string,gen,11)==0)
					{
						sscanf(string,"..generated %d conformers",&conf_number);
						if(total_conf!=0)
						printf("round %d total conformations %d generation %d conformers ratio %.3f%%\n",round,total_conf,conf_number,((conf_number*100.0)/total_conf));
						if((total_conf==0)&&(rotbond==0)){
							printf("NO conformation generated, End conformational search");
							fclose(log_fp);
							remove(out_name);
							remove("confab.log");
							return 1;
						}
						if(total_conf<=500){
						if( conf_number==total_conf  )
						{	
							fclose(log_fp);
							return 0;
						}
						else
						{
							fclose(log_fp);
							remove(out_name);
							remove("confab.log");
						}							
							
						}else{
							if( (conf_number>=500)||    (((conf_number*100.0)/total_conf)>=10.0)    )
							{	
								fclose(log_fp);
								return 0;
							}
							else
							{
								fclose(log_fp);
								remove(out_name);
								remove("confab.log");
							}							
							
							
						}

					}

				}
			}
			else{
					printf("\n Open confab.log failed\n") ;
					return -1;
				}
		
		
	}


}


int cal_conformer(char *file_name)
{
	
	//Obabel  9.mol2 -O        9_4.xyz  --conformer  --nconf 200  --writeconformers   --random  
	if(is_obabel_installed()!=1){
		printf("The presence of obabel was not detected on this device.\n Please install, recommended version: Open Babel 3.0.0\n");
		exit( -1);	
	}
	char out_name[50],call_confab[200];
	
	out_name[0]='\0';
	strcat(out_name,"cos_");
	strcat(out_name,file_name);
	int i=0;
	while(out_name[i]!='\0') i++;
	while(out_name[i]!='.') i--;
	out_name[i]='\0';
	strcat(out_name,".xyz");	
	
	remove(out_name);

	call_confab[0]='\0';
	strcat(call_confab,"obabel ");
	strcat(call_confab,file_name);
	strcat(call_confab,"  -O  ");
	strcat(call_confab,out_name);
	strcat(call_confab," --conformer  --nconf 500  --writeconformers   --random  ");
	printf("%s\n",call_confab);
	system(call_confab);
	return 0;		
}

int cal_obabel_mmf(char *file_name)
{
	if(file_type(file_name)!=4){
		printf("\n\n\n    ERROR The input file [%s] for energy minimization should be in xyz format\n\n",file_name);
		exit(-1);
	}

	int atom_number=0,file_number=0;
	int file_count=0;
	char mm_name[100],command[300];
	mm_name[0]='\0';
	
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);



//----------out result-----------
	strcat(mm_name,"mm_");
	strcat(mm_name,file_name);	
	int i=0;
	while(mm_name[i]!='\0') i++;
	while(mm_name[i]!='.') i--;
	mm_name[i]='\0';
	strcat(mm_name,".xyz");
//--------------------------------
	if(file_number==1){
		// Obabel  1.xyz   -O       1_opt.xyz    --minimize  --append "Energy"  --ff MMFF94   --steps 2500
		command[0]='\0';
		strcat(command,"Obabel ");
		strcat(command,file_name);
		strcat(command," -O ");
		strcat(command,mm_name);
		strcat(command," --minimize  --append \"Energy\"  --ff MMFF94   --steps 3000");
		printf("%s",command);
		system(command);			
	}else if(file_number>1){
	//	remove(mm_name);
		FILE *mmresult;
		mmresult=fopen(mm_name, "w");
		if(mmresult==NULL){
			printf("Open file %s fialed\n",mm_name);
			exit(-1);	
		}
//-----------------------METHOD2  trace file-----------------------------------------

	FILE *file;
	char fs[300];
	file=fopen(file_name, "r");
	if(file==NULL){
		printf("Open file %s fialed\n",file_name);
		exit(-1);	
	}
	int cf=0;
	FILE *temp,*mm_temp;
	temp=fopen("temp.xyz","w");
	if(temp==NULL){
		printf("Open file temp.xyz fialed\n");
		fclose(file);
		fclose(mmresult);
		exit(-1);
	}
	
	char fstring[300];	
	while(fgets(fs,255,file)!=NULL){

		
		
		cf++;
		
		
		
		if((cf-((atom_number+2)*file_count))==2){
			fputs("\n",temp);	
		}
		else{
			fputs(fs,temp);
		}
		
		if((cf%(atom_number+2))==0){
			file_count++;
			fclose(temp);
			
			command[0]='\0';
			strcat(command,"Obabel ");
			strcat(command,"temp.xyz");
			strcat(command," -O ");
			strcat(command,"mm_temp.xyz");
			strcat(command," --minimize  --append \"Energy\"  --ff MMFF94   --steps 3000 ");
			//Hide command output and display progress bar
			#ifdef _WIN32
				strcat(command," > nul 2>&1");
			#else
				strcat(command," > /dev/null 2>&1");
			#endif			
			// Display progress bar 
			print_progress_bar1(file_count, file_number, "MMFF94 energy minimization");			
			
			
		//	printf("Frame %d of %d ",file_count,file_number);
		//	printf("%s",command);
			system(command);
			FILE *mm_temp;
			mm_temp=fopen("mm_temp.xyz","r");
			if(mm_temp!=NULL){
				while(fgets(fstring,255,mm_temp)!=NULL){
					fputs(fstring, mmresult);
				}
				fclose(mm_temp);
			}
			else{
				printf("Open file mm_temp.xyz fialed\n");
				fclose(mmresult);
				fclose(file);
				exit(-1);
			}
			if(file_count<=(file_number-1))
			temp=fopen("temp.xyz","w");
			if(temp==NULL){
				printf("Open file temp.xyz fialed\n");
				fclose(mmresult);
				fclose(file);
				exit(-1);
			}		
	
		}
		
	}

//-------------------------------------------------------------
		fclose(mmresult);
		remove("mm_temp.xyz");
		remove("temp.xyz");
		
		
	
		
	}
	
	return 0;	
}





int get_mopac_dataf(char *file_name,char *out_name)
{
	char *energy_flag="          FINAL HEAT OF FORMATION =",
		 fstring[200],
		 *zuobiao="                             CARTESIAN COORDINATES",
		 *atom_number="           Empirical Formula: ",
		 atom_arry[5], *mopac_done=" == MOPAC DONE ==";
	FILE *fp,*out_fp;
	int count=0,zuobiao_flag=0,blank_line=0,i=0,atom=0,finished=0;
	float energy_kcal=0.0,energy_kj=0.0,x,y,z;
/*******************************if calculate finished**********************************/	
	fp=fopen(file_name,"r");
		if(fp!=NULL)
	{
		while ((fgets(fstring,199,fp))!=NULL)
			{

				if(strncmp(fstring,mopac_done,17)==0)
				{
					finished=1;
					break;
				}
				
								
			}
	}
	else
	{
		printf ("\n Open file %s failed \n",file_name);
		exit(-1);		
	}
	fclose(fp);	
	if(finished==0) return -1;
	
	
	
	
/*******************************out xyz name**********************************/


	out_fp=fopen(out_name,"w");

/*****************************************************************/	
	
	
	
	
	
	
	
/**********************get atom number ****************************/
	fp=fopen(file_name,"r");
		if(fp!=NULL)
	{
		while ((fgets(fstring,199,fp))!=NULL)
			{
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				if(strncmp(fstring,atom_number,29)==0)
				{
					i=0;
					while(fstring[i]!='=')
					{
						fstring[i]=' ';
						i++;
					}
					sscanf(fstring,"                                            = %d atoms",&atom);
					break;
				}
				
								
			}
	}
	else
	{
		printf ("\n Open file %s failed \n",file_name);
		exit(-1);		
	}
	

fprintf(out_fp,"%d\n",atom);
/*************************************************/	
	rewind(fp);
	if(fp!=NULL)
	{
		while ((fgets(fstring,199,fp))!=NULL)
			{
				if (strlen(fstring) > 1 && (fstring[strlen(fstring)-2]=='\r')){
      				fstring[strlen(fstring)-2]='\n';
      				fstring[strlen(fstring)-1]='\0';
    			}
				
				if(strncmp(fstring,energy_flag,35)==0)
				{
					sscanf(fstring,"          FINAL HEAT OF FORMATION =        %f KCAL/MOL =    %f KJ/MOL",&energy_kcal,&energy_kj);
					fprintf(out_fp,"%s %-15.7f\n","Energy(KJ/MOL)",energy_kj);
				//	printf("%.5f KJ/MOL %.5f KCAL/MOL\n",,energy_kcal);
				}
				if(strncmp(fstring,zuobiao,35)==0)
				{
					zuobiao_flag=1;
				}
				if((zuobiao_flag==1)&&(strncmp(fstring,zuobiao,35)!=0))
				{
					if(strlen(fstring)==1) blank_line++;
					if((strlen(fstring)!=1)&&(blank_line<=1))
					{
						i=0;
						while(isalpha(fstring[i])==0)
						{
							fstring[i]=' ';
							i++;
						}
						sscanf(fstring,"        %s       %f    %f    %f",atom_arry,&x,&y,&z);
						//printf("[a]%s",fstring);
    					fprintf(out_fp,"%3s%15.5f%15.5f%15.5f\n",atom_arry,x,y,z);
  																	
					}

				}								
			}
	}
	else
	{
		printf ("\n Open file %s failed \n",file_name);
		exit(-1);		
	}
	fclose(fp);	
	fclose(out_fp);
	return 0;
}




void print_progress_bar(int progress, int total) {
    int barWidth = 70; 
    float progressPercent = (float)progress / total;
    int pos = barWidth * progressPercent,i;
 
    printf("[");
    for (i = 0; i < barWidth; ++i) {
        if (i < pos) {
            printf("=");
        } else if (i == pos) {
            printf(">");
        } else {
            printf(" ");
        }
    }
    printf("] %d %%\r", (int)(progressPercent * 100));
    fflush(stdout);
}






 
int cal_mopac_mmf(char *file_name,int multiplicity, int charge, float permittivity,char *mopac_path )
{
	if(file_type(file_name)!=4){
		printf("\n\n\n    ERROR The input file for energy minimization should be in xyz format\n\n");
		exit(-1);
	}	

	int atom_number=0,file_number=0;
	
	char mm_name[100],command[300];
	mm_name[0]='\0';
	

	
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);

	
	
	
		
	if(file_number==1){
		
		mm_name[0]='\0';
		strcat(mm_name,"mopac_");
		strcat(mm_name,file_name); 
		int i=0;
		while(mm_name[i]!='\0') i++;
		while(mm_name[i]!='.') i--;
		mm_name[i]='\0';
		strcat(mm_name,".xyz");
		FILE *mopac_temp; 
		mopac_temp=fopen("temp_mopac.mop","w");
		if(mopac_temp==NULL){
			printf("Open file mopac_temp.mop fialed\n");
			exit(-1);			
		}
		else{
			fprintf(mopac_temp,"CHARGE=%d ",charge);
			
			if(multiplicity==1) fprintf(mopac_temp," SINGLET ");
			else if(multiplicity==2) fprintf(mopac_temp," DOUBLET ");
			else if(multiplicity==3) fprintf(mopac_temp," TRIPLET ");
			else if(multiplicity==4) fprintf(mopac_temp," QUARTET ");
			else if(multiplicity==5) fprintf(mopac_temp," QUINTET ");
			else{
				printf("Incorrect multiplicity parameter, please carefully check.\n This method only supports: 1-5.\n");
				fclose(mopac_temp);
				exit(-1);
			}
			
			fprintf(mopac_temp,"PM6-D3H4 precise OPT");
			if((permittivity>1.0)&&(permittivity<200.0))
			fprintf(mopac_temp," eps=%f ",permittivity);
			fprintf(mopac_temp,"\n");
			fprintf(mopac_temp,"molecule\n");
			fprintf(mopac_temp,"All coordinates are Cartesian\n");	
			int i;
			
			FILE *input_filefp; 
			input_filefp=fopen(file_name,"r");
			if(input_filefp==NULL){
				printf("Open file %s fialed\n",file_name);
				exit(-1);			
			}else{
				char fs[300], ele[5];
				int count=0;
				float x,y,z;
				while ((fgets(fs,256,input_filefp))!=NULL){
					if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      					fs[strlen(fs)-2]='\n';
      					fs[strlen(fs)-1]='\0';
    				}
					count++;
					if(count>2){
						sscanf(fs,"%s %f %f %f",ele,&x,&y,&z);
						fprintf(mopac_temp,"%-2s %10.6f  1  %10.6f  1  %10.6f  1  \n",ele,x,y,z);
					}
					
				}
				fclose(input_filefp);
			}			
			fclose(mopac_temp);
			command[0]='\0';
			strcat(command,mopac_path);
			strcat(command," ");
			strcat(command,"temp_mopac.mop");
		//	printf("%s",command);
			system(command);
			if(get_mopac_dataf("temp_mopac.out",mm_name)!=0){
				printf("Failed to run mopac\n");
				exit(-1);	
			}	
		}
		
	}else if(file_number>1){
		
		int i=0,j=0;
		mm_name[0]='\0';
		strcat(mm_name,"mopac_");
		strcat(mm_name,file_name); 
		while(mm_name[i]!='\0') i++;
		while(mm_name[i]!='.') i--;
		mm_name[i]='\0';
		strcat(mm_name,".xyz");		
		FILE *mmresult;
		mmresult=fopen(mm_name,"w");

		if(mmresult==NULL){
			printf("Open file fialed %s\n",mm_name);
			exit(-1);
		}
		
		FILE *mopac_temp; 
		mopac_temp=fopen("temp_mopac.mop","w");
		if(mopac_temp==NULL){
			printf("Open file mopac_temp.mop fialed\n");
			exit(-1);			
		}else{
			fprintf(mopac_temp,"CHARGE=%d ",charge);
			
			if(multiplicity==1) fprintf(mopac_temp," SINGLET ");
			else if(multiplicity==2) fprintf(mopac_temp," DOUBLET ");
			else if(multiplicity==3) fprintf(mopac_temp," TRIPLET ");
			else if(multiplicity==4) fprintf(mopac_temp," QUARTET ");
			else if(multiplicity==5) fprintf(mopac_temp," QUINTET ");
			else{
				printf("Incorrect multiplicity parameter, please carefully check.\n This method only supports: 1-5.\n");
				fclose(mopac_temp);
				exit(-1);
			}
			
			
			
			fprintf(mopac_temp,"PM6-D3H4 precise OPT");
			if((permittivity>1.0)&&(permittivity<200.0))
			fprintf(mopac_temp," eps=%f ",permittivity);
			fprintf(mopac_temp,"\n");
			fprintf(mopac_temp,"molecule\n");
			fprintf(mopac_temp,"All coordinates are Cartesian\n");		
		}
		
		
		FILE *file;
		char fs[300];
		file=fopen(file_name, "r");
		if(file==NULL){
			printf("Open file %s fialed\n",file_name);
			exit(-1);	
		}
		int line=0,file_count=0; 
		float x,y,z;
		char ele[5];
		while(fgets(fs,255,file)!=NULL){
			if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      			fs[strlen(fs)-2]='\n';
      			fs[strlen(fs)-1]='\0';
    		}
			
			line++;
			if((line-((atom_number+2)*file_count))>2){
				sscanf(fs,"%s %f %f %f",ele,&x,&y,&z);
				fprintf(mopac_temp,"%-2s %10.6f  1  %10.6f  1  %10.6f  1  \n",ele,x,y,z);					
			}
			if((line%(atom_number+2))==0){
				file_count++;
				fclose(mopac_temp);
				command[0]='\0';
				
				strcat(command,mopac_path);
				strcat(command," ");
				strcat(command,"temp_mopac.mop");
				printf("%s",command);
				system(command);				
				printf("Running progress:\n");
				print_progress_bar((file_count*100/file_number), 100);
				
				if(get_mopac_dataf("temp_mopac.out","temp_mopac.xyz")!=0){
					printf("Failed to run mopac\n");
					fclose(mmresult);
					exit(-1);	
				}
				FILE *temp_mopac_result;
				temp_mopac_result=fopen("temp_mopac.xyz","r");
				if(temp_mopac_result==NULL){
					printf("Open file fialed temp_mopac.xyz\n");
					fclose(mmresult);
					exit(-1);				
				}else{
					char fs[300];
					while(fgets(fs,255,temp_mopac_result)!=NULL){
						fputs(fs, mmresult);
					}
					fclose(temp_mopac_result);	
					remove("temp_mopac.xyz");
				}				
				
				if(file_count<=(file_number-1)){
					mopac_temp=fopen("temp_mopac.mop","w");
					if(mopac_temp==NULL){
						printf("Open file mopac_temp.mop fialed\n");
						exit(-1);			
					}else{
						fprintf(mopac_temp,"CHARGE=%d ",charge);
						fprintf(mopac_temp,"PM6-D3H4 precise OPT");
						if((permittivity>1.0)&&(permittivity<200.0))
						fprintf(mopac_temp," eps=%f ",permittivity);
						fprintf(mopac_temp,"\n");
						fprintf(mopac_temp,"molecule\n");
						fprintf(mopac_temp,"All coordinates are Cartesian\n");		
					}	
				}

			}

			
		}
		fclose(mmresult);
	
	}else{
		printf("ERROR read file number failed %d\n",file_number);
	}
	
	printf("\n\n %s, file MOPAC energy minimization completed\n",file_name);
	return 0;
}


int get_energy(char *file_name,float bloz)
{
	float blozsel;
	if(bloz==0.0) blozsel=0.5;
	else blozsel=bloz;
	int atom_number,file_number,count=0,fcount=0,linenumber=0;
	FILE *fp;
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);	
	if((atom_number<=1)||(file_number<=1)){
		printf("\n\n read arom number %d or file number %d failed of file %s\n",atom_number,file_number,file_name);
		return 1;
	//	exit(-1);	
	}
	struct energy{
		int id;
		float energy;
	};
	struct energy  *energy_idx;
	energy_idx=(struct energy *)malloc(file_number*sizeof(struct energy));
	char fs[300];
	fp=fopen(file_name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",file_name);
		return -1;
	}else{
		while(fgets(fs,255,fp)!=NULL){
			if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      			fs[strlen(fs)-2]='\n';
      			fs[strlen(fs)-1]='\0';
    		}
			
			
 			linenumber++;
 			if((linenumber-((atom_number+2)*fcount))==2){
 				float energy;
 				int i=0;
 				while(fs[i]!=' ') {
 					fs[i]=' ';
 					i++;
				 }
 				sscanf(fs,"%f",&energy);
 				energy_idx[fcount].energy=energy;
 				energy_idx[fcount].id=fcount+1;
			 }
			 if((linenumber%(atom_number+2))==0){
			 	fcount++;
			 }
		}	
		fclose(fp);
	}
	int tempid;
	float tempe;
	int i,j;
	for(i=0;i<=file_number-1;i++){
		for(j=i+1;j<=file_number-1;j++){
			if(energy_idx[i].energy>=energy_idx[j].energy){
				tempid=energy_idx[i].id;
				tempe=energy_idx[i].energy;
				energy_idx[i].id=energy_idx[j].id;
				energy_idx[i].energy=energy_idx[j].energy;
				energy_idx[j].id=tempid;
				energy_idx[j].energy=tempe;
			}
		}
	}
	float tempature=298.15,R=0.0019872,dis=0;
	
		for(i=0;i<=file_number-1;i++)
	{
		dis=exp(-(energy_idx[i].energy-energy_idx[0].energy)/(tempature*R))+dis;
	}
	
	printf("\n\nframe  Energy  Distribution(%%)\n");
	
	// choose boltzmann distribution >0.5%
	for(i=0;i<=file_number-1;i++)
	{
		printf("%-5d  %-6.3f  %-6f%%\n",energy_idx[i].id,energy_idx[i].energy,exp(-(energy_idx[i].energy-energy_idx[0].energy)/(tempature*R))/dis*100);
		if((exp(-(energy_idx[i].energy-energy_idx[0].energy)/(tempature*R))/dis*100)<blozsel){
		energy_idx[i].id=-1;	
		}
		
	}
	
	char choosed[100];
	choosed[0]='\0';
	strcat(choosed,file_name);
	i=0;
	while(choosed[i]!='\0') i++;
	while(choosed[i]!='.') i--;
	choosed[i]='\0';
	strcat(choosed,"_choosed.xyz");
	FILE *choose_fp;
	choose_fp=fopen(choosed,"w");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",choosed);
		return -1;
	}
	
	fp=fopen(file_name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",file_name);
		return -1;
	}else{
		
		for(i=0;i<=file_number-1;i++){
			rewind(fp);
			if(energy_idx[i].id>=1){
				linenumber=0;
				while(fgets(fs,255,fp)!=NULL){
 					linenumber++;
 					if((linenumber>((atom_number+2)*(energy_idx[i].id-1)))&&(linenumber<=((atom_number+2)*(energy_idx[i].id)))){
 						fputs(fs,choose_fp);
					 }
				}	
			}
		}
		fclose(fp);
		fclose(choose_fp);
		
	}	
	
	
	
	
	free(energy_idx);
	return 0;
}
 

int get_dist_data(char *file_name,float *distribution, int input_filenumber )
{

	int atom_number,file_number,count=0,fcount=0,linenumber=0;
	float sel_bloz;
	FILE *fp;
	atom_number=get_xyz_atom_number(file_name);
	file_number=get_xyz_file_number(file_name);	
	if(input_filenumber!=file_number){
		printf("\n\nThe number of files entered %d does not match the number %d of files in %s file",input_filenumber,file_number,file_name);
		exit(-1);
	}
	if((atom_number<=1)||(file_number<=1)){
		printf("\n\n read atom number %d or file number %d failed of file %s\n",atom_number,file_number,file_name);
		exit(-1);	
	}
	struct energy{
		int id;
		float energy;
	};
	struct energy  *energy_idx;
	energy_idx=(struct energy *)malloc(file_number*sizeof(struct energy));
	char fs[300];
	fp=fopen(file_name,"r");
	if(fp==NULL){
		printf("\n\n %s, file failed\n",file_name);
		return -1;
	}else{
		while(fgets(fs,255,fp)!=NULL){
			if (strlen(fs) > 1 && (fs[strlen(fs)-2]=='\r')){
      			fs[strlen(fs)-2]='\n';
      			fs[strlen(fs)-1]='\0';
    		}
			
			
 			linenumber++;
 			if((linenumber-((atom_number+2)*fcount))==2){
 				float energy;
 				int i=0;
 				while(fs[i]!=' ') {
 					fs[i]=' ';
 					i++;
				 }
 				sscanf(fs,"%f",&energy);
 				energy_idx[fcount].energy=energy;
 				energy_idx[fcount].id=fcount+1;
			 }
			 if((linenumber%(atom_number+2))==0){
			 	fcount++;
			 }
		}	
		fclose(fp);
	}
	int tempid;
	float min_energy=energy_idx[0].energy;
	int i;
	for(i=0;i<=file_number-1;i++){
		if(min_energy>=energy_idx[i].energy)
		min_energy=energy_idx[i].energy; 

	}
	float tempature=298.15,R=0.0019872,dis=0;
	
	for(i=0;i<=file_number-1;i++)
	{
		dis=exp(-(energy_idx[i].energy-min_energy)/(tempature*R))+dis;
	}
	
	for(i=0;i<=file_number-1;i++)
	{
		distribution[i]=exp(-(energy_idx[i].energy-min_energy)/(tempature*R))/dis*100;
	}
	free(energy_idx);
	return 0;
}
  


int ask_if_bloz()
{
	int blot;
	char buffer_command[300];
	printf("\nDo we need to perform Boltzmann distribution screening first:enter ""1 (Yes) "" "" 0 (No)""\n" );
	printf("Please enter the selected method:\n");
	fgets(buffer_command, sizeof(buffer_command), stdin);
	sscanf(buffer_command, "%d", &blot);
	return blot;
}
float get_dis()
{
	char buffer_command[300];
	float bloz;
	printf("Enter the lowest proportion of the selected Boltzmann distribution 0-100 (0.5 refers to 0.5%% and default):\n");
	fgets(buffer_command, sizeof(buffer_command), stdin);
	sscanf(buffer_command, "%f", &bloz);
	printf("You selected Boltzmann distribution is：%f%%\n",bloz);
	return bloz;	
}
void add_prefix(char *newnamef, char *name, char *add)
{
	int i=0;
	newnamef[0]='\0';
	strcat(newnamef,add);
	strcat(newnamef,name);
	while(newnamef[i]!='\0') i++;
	while(newnamef[i]!='.') i--;	
	return ;
}
void add_suffix(char *newnamef, char *name, char *add)
{
	int i=0;
	newnamef[0]='\0';
	strcat(newnamef,name);
	while(newnamef[i]!='\0') i++;
	while(newnamef[i]!='.') i--;
	newnamef[i]='\0';
	strcat(newnamef,add);
	return ;
}



void set_name(int cal_type,char *input_name,char *seted_name)
{
//	printf("calculate type is %d\n\n",cal_type);
	char input_file[100],newname[100], newname1[100], newname2[100], newname3[100], newname4[100], newname5[100];
	char name_cal[100];
	int i=0;
	input_file[0]='\0';
	if(cal_type!=1){
		strcat(input_file,input_name);
		i=0;
		while(input_file[i]!='\0'){
			i++;
		}
		while(input_file[i]!='.'){
			i--;
		}	
		input_file[i]='\0';
		strcat(input_file,".xyz");			
	}

	name_cal[0]='\0';
	if(cal_type==1){
		strcat(name_cal,input_name);
	}else if(cal_type==2){  // only conformational search
		add_prefix(name_cal, input_file, "cos_");
	}else if(cal_type==3){    // only MMFF94 bloz
		add_prefix(newname, input_file, "mm_");
		add_suffix(name_cal, newname, "_choosed.xyz");
	}else if(cal_type==4){          // only MMFF94 
		add_prefix(name_cal, input_file, "mm_");
	}else if(cal_type==5){ // only mopac bloz
		add_prefix(newname, input_file, "mopac_");
		add_suffix(name_cal, newname, "_choosed.xyz");		
	}else if(cal_type==6){ // only mopac
		add_prefix(name_cal, input_file, "mopac_");		
	}else if(cal_type==7){ // mmf94 bloz + mopac bloz	
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz		
		add_suffix(newname1, newname, "_choosed.xyz"); // bloz mm_a_choose.xyz
		add_prefix(newname2, newname1, "mopac_");  // mopac mopac_mm_a_choose.xyz
		add_suffix(name_cal, newname2, "_choosed.xyz");	// bloz mopac_mm_a_choose_choose.xyz	
	}else if(cal_type==8){
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz		
		add_suffix(newname1, newname, "_choosed.xyz"); // bloz mm_a_choose.xyz
		add_prefix(name_cal, newname1, "mopac_");  // mopac mopac_mm_a_choose.xyz		
	}else if(cal_type==9){
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz	
		add_prefix(newname1, newname, "mopac_");  // mopac mopac_mm_a.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz mm_a_choose.xyz		
	}else if(cal_type==10){
		add_prefix(newname, input_file, "mm_");  // mmff94 mm_a.xyz	
		add_prefix(name_cal, newname, "mopac_");  // mopac mopac_mm_a.xyz
	}else if(cal_type==11){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz mm_a_choose.xyz		
	}else if(cal_type==12){
		add_prefix(newname, input_file, "cos_");  // cof cos_a.xyz	
		add_prefix(name_cal, newname, "mm_");  // mmff94 mm_cos_a.xyz
	}else if(cal_type==13){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mopac_");  // mopac mopac_cos_a.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz mopac_cos_a_choose.xyz		
	}else if(cal_type==14){
		add_prefix(newname, input_file, "cos_");  // cof cos_a.xyz	
		add_prefix(name_cal, newname, "mopac_");  // mopac mopac_cos_a.xyz
	}else if(cal_type==15){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_suffix(newname2, newname1, "_choosed.xyz"); // bloz mm_cos_a_choose.xyz	
		add_prefix(newname3, newname2, "mopac_");  // mopac mopac_mm_cos_a_choose.xyz
		add_suffix(name_cal, newname3, "_choosed.xyz"); // bloz mopac_mm_cos_a_choose_choose.xyz		
	}else if(cal_type==16){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_suffix(newname2, newname1, "_choosed.xyz"); // bloz mm_cos_a_choose.xyz	
		add_prefix(name_cal, newname2, "mopac_");  // mopac mopac_mm_cos_a_choose.xyz
	}else if(cal_type==17){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_prefix(newname2, newname1, "mopac_");  // mopac mopac_mm_cos_a.xyz
		add_suffix(name_cal, newname2, "_choosed.xyz"); // bloz mopac_mm_cos_a_choose.xyz		
	}else if(cal_type==18){
		add_prefix(newname, input_file, "cos_");  // cnf cos_a.xyz	
		add_prefix(newname1, newname, "mm_");  // mmff94 mm_cos_a.xyz
		add_prefix(name_cal, newname1, "mopac_");  // mopac mopac_mm_cos_a.xyz
	}else if(cal_type==19){
		add_suffix(newname, input_file, "_choosed.xyz"); // bloz a_choose.xyz
		add_prefix(newname1, newname, "mopac_");  // mopac mopac_a_choose.xyz
		add_suffix(name_cal, newname1, "_choosed.xyz"); // bloz a_choose_choose.xyz
	}else if(cal_type==20){
		add_prefix(newname, input_file, "mopac_");
		add_suffix(name_cal, newname, "_choosed.xyz");
	}else{
		printf("Error name to calculate\n");
		exit(-1);
	}
	seted_name[0]='\0';
	strcat(seted_name,name_cal);
	
	
	return ;
	
}

void chang_suffix2xyz(char *name)
{
	int i=0;
	while(name[i]!='\0') i++;
	while(name[i]!='.') i--;
	name[i]='\0';
	strcat(name,".xyz");
	return ;
}


// The entire string can be parsed as a number
int is_number(const char *str) {
    char *endptr;
    strtod(str, &endptr);
    return (*endptr == '\0');  
}







void cal_3d_hu_moment(struct hu_data *hdata, struct grid  *grid_data, int grid_count)
{	
	int i;
	double m000=0.00,m100=0.00,m010=0.00,m001=0.00,x_c=0.00,y_c=0.00,z_c=0.00,u200=0.00,u020=0.00,u002=0.00,ks=0.00,u000=0.00;
	double u110=0.00,u101=0.00,u011=0.0,u300=0.00,u030=0.00,u003=0.00;
	double u120=0.00,u102=0.00,u012=0.00,u210=0.00,u021=0.00,u201=0.00,u111=0.0;
	double tt1=0.0,tt2=0.0,tt3=0.0,tt4=0.0;
	for(i=0;i<=grid_count-1;i++)
	{	
		if  (grid_data[i].state==1)
		{    
			m000=m000+1.000;
			m100=m100+grid_data[i].x;
			m010=m010+grid_data[i].y;
			m001=m001+grid_data[i].z;
		}	
	}	
	x_c=m100/m000;
	y_c=m010/m000;
	z_c=m001/m000;

	for(i=0;i<=grid_count-1;i++)
	{	
		if  (grid_data[i].state==1)
		{    
			u200=u200+((grid_data[i].x-x_c)*(grid_data[i].x-x_c));
			u020=u020+((grid_data[i].y-y_c)*(grid_data[i].y-y_c));
			u002=u002+((grid_data[i].z-z_c)*(grid_data[i].z-z_c));
			u110=u110+((grid_data[i].x-x_c)*(grid_data[i].y-y_c));
			u101=u101+((grid_data[i].x-x_c)*(grid_data[i].z-z_c));
			u011=u011+((grid_data[i].y-y_c)*(grid_data[i].z-z_c));
			u300=u300+((grid_data[i].x-x_c)*(grid_data[i].x-x_c)*(grid_data[i].x-x_c));
			u030=u030+((grid_data[i].y-y_c)*(grid_data[i].y-y_c)*(grid_data[i].y-y_c));
			u003=u003+((grid_data[i].z-z_c)*(grid_data[i].z-z_c)*(grid_data[i].z-z_c));			
			u120=u120+((grid_data[i].x-x_c)*(grid_data[i].y-y_c)*(grid_data[i].y-y_c));
			u102=u102+((grid_data[i].x-x_c)*(grid_data[i].z-z_c)*(grid_data[i].z-z_c));
			u012=u012+((grid_data[i].y-y_c)*(grid_data[i].z-z_c)*(grid_data[i].z-z_c));
			u210=u210+((grid_data[i].x-x_c)*(grid_data[i].x-x_c)*(grid_data[i].y-y_c));
			u021=u021+((grid_data[i].y-y_c)*(grid_data[i].y-y_c)*(grid_data[i].z-z_c));
			u201=u201+((grid_data[i].x-x_c)*(grid_data[i].x-x_c)*(grid_data[i].z-z_c));
			u111=u111+((grid_data[i].x-x_c)*(grid_data[i].y-y_c)*(grid_data[i].z-z_c));
			
		}	
	}
	u000=m000;

			u200=u200/(pow(u000,(5.00000/3.0000)));
			u020=u020/(pow(u000,(5.00000/3.0000)));
			u002=u002/(pow(u000,(5.00000/3.0000)));
			u110=u110/(pow(u000,(5.00000/3.0000)));
			u101=u101/(pow(u000,(5.00000/3.0000)));
			u011=u011/(pow(u000,(5.00000/3.0000)));
			u300=u300/(pow(u000,2.00000));
			u030=u030/(pow(u000,2.0000));
			u003=u003/(pow(u000,2.0000));		
			u120=u120/(pow(u000,2.0000));
			u102=u102/(pow(u000,2.0000));
			u012=u012/(pow(u000,2.0000));
			u210=u210/(pow(u000,2.0000));
			u021=u021/(pow(u000,2.0000));
			u201=u201/(pow(u000,2.0000));
			u111=u111/(pow(u000,2.0000));	
	
//	
	
	
	
//	t1s=(u200+u020+u002)/(pow(u000,(5.00000/3.0000)));

	
	
	
	tt1=u200+u020+u002;
	tt2=(u200*u020)+(u200*u002)+(u020*u002)-(u110*u110)-(u101*u101)-(u011*u011);
	tt3=(u200*u020*u002)+(2.0*u110*u101*u011)-(u200*u011*u011)-(u020*u101*u101)-(u002*u110*u110);
	tt4=(u300*u300)+(u030*u030)+(u003*u003)+(3.0*u120*u120)+(3.0*u102*u102)+(3.0*u012*u012)+(3.0*u210*u210)+(3.0*u021*u021)+(3.0*u201*u201)+(6.0*u111*u111);
//	printf(" h1 %f h2 %f h3 %f h4 %f",tt1,tt2,tt3,tt4);
	hdata->v4=tt1;
	hdata->v5=tt2;
	hdata->v6=tt3;
	hdata->v7=tt4;
	
//	ks =(pow(3.0,(5.0/3.0))/(5.0*(pow((4.0*3.1415926),(2.0/3.0)))))*(1.0/tt1);
	ks = pow(3.0, 5.0/3.0) / (5.0 * pow(4.0 * 3.1415926, 2.0/3.0)) * (1.0/tt1);

//	printf(" ks %f",ks);
	hdata->v8=ks;	
	return ;
	
	
}




















int cal_moment(struct hu_data *data, struct grid  *grid_data, int grid_count,float grid_gap, float volume, float suface)
{

	data->v3=suface/volume;
	cal_3d_hu_moment(data, grid_data, grid_count);
	double kst=0.00;
	kst=(36.000000*3.14159265358*volume*volume)/(suface*suface*suface);
	data->v9=kst;
	double vr=0.0,vinsph=0.0,kfit=0.0;
	long int dot_in_sphere=0,dot_in_mol=0;
	float cx=0,cy=0,cz=0;
	int i;
	for(i=0;i<=grid_count-1;i++)
	{
		if((grid_data[i].state==2)|| (grid_data[i].state==1))
		{
			dot_in_mol++;
			cx=cx+grid_data[i].x;
			cy=cy+grid_data[i].y;
			cz=cz+grid_data[i].z;
		}	
	}
	cx=cx/dot_in_mol;
	cy=cy/dot_in_mol;
	cz=cz/dot_in_mol;
	vr=pow((volume*3.0/(4.0*3.14159265358)),(1.0/3.0));
	for(i=0;i<=grid_count-1;i++)
	{
		if((grid_data[i].state==2)|| (grid_data[i].state==1))
		{
			if(vr>=(sqrt(  (grid_data[i].x-cx)*(grid_data[i].x-cx) + (grid_data[i].y-cy)*(grid_data[i].y-cy)  + (grid_data[i].z-cz)*(grid_data[i].z-cz))))
			{
				dot_in_sphere++;
			}
		}	
	}
	vinsph=dot_in_sphere*grid_gap*grid_gap*grid_gap;
	kfit=vinsph/((2.0*volume)-vinsph);
	data->v10=kfit;
	return 0;
}







//Extract PDB ID from composite name 
void extract_pdb_id(const char *compound_name, char *pdb_id) {
    //Check if it starts with "L_" and has sufficient length 
    if (strncmp(compound_name, "L_", 2) == 0 && strlen(compound_name) >= 6) {
        //Copy the 3rd to 6th characters directly (index 2-5) 
        strncpy(pdb_id, compound_name + 2, 4);
        pdb_id[4] = '\0';
        
        //Ensure it is capitalized 
        int i;
        for ( i = 0; pdb_id[i] != '\0'; i++) {
            pdb_id[i] = toupper(pdb_id[i]);
        }
    } else {
        //If it is not in standard format, return an empty string
        pdb_id[0] = '\0';
    }
}


 






   






//Read the Targets. tsv file
int read_targets_info(const char *filename, struct target_info **targets_ptr) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Could not open file %s\n", filename);
        *targets_ptr = NULL;
        return 0;
    }
    
    //Calculate the number of rows for the first time
    int count = 0;
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        count++;
    }
    
    //Allocate memory (subtract header row) 
    *targets_ptr = (struct target_info*)malloc((count-1) * sizeof(struct target_info));
    if (*targets_ptr == NULL) {
        printf("Memory allocation failed for targets info\n");
        fclose(file);
        return 0;
    }
    
    //Second reading of data 
    rewind(file);
    int actual_count = 0;
    
    //Skip title line 
    fgets(line, sizeof(line), file);
    
    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;
        
        //Skip blank lines
        if (strlen(line) == 0) continue;
        
        //Split strings using strtok (tab separated) 
        char *token = strtok(line, "\t");
        if (token == NULL) continue;
        
        //Copy PDB_ID 
        strncpy((*targets_ptr)[actual_count].pdb_id, token, 
                sizeof((*targets_ptr)[actual_count].pdb_id) - 1);
        (*targets_ptr)[actual_count].pdb_id[sizeof((*targets_ptr)[actual_count].pdb_id) - 1] = '\0';
        
        //Read protein names
        token = strtok(NULL, "\t");
        if (token == NULL) continue;
        strncpy((*targets_ptr)[actual_count].protein_names, token, 
                sizeof((*targets_ptr)[actual_count].protein_names) - 1);
        (*targets_ptr)[actual_count].protein_names[sizeof((*targets_ptr)[actual_count].protein_names) - 1] = '\0';
        
        //Read Gene Names 
        token = strtok(NULL, "\t");
        if (token == NULL) continue;
        strncpy((*targets_ptr)[actual_count].gene_names, token, 
                sizeof((*targets_ptr)[actual_count].gene_names) - 1);
        (*targets_ptr)[actual_count].gene_names[sizeof((*targets_ptr)[actual_count].gene_names) - 1] = '\0';
        
        //Read Uniprot ID (new field) 
        token = strtok(NULL, "\t");
        if (token == NULL) {
            //If Uniprot ID is missing, set to an empty string 
            (*targets_ptr)[actual_count].uniprot_id[0] = '\0';
        } else {
            strncpy((*targets_ptr)[actual_count].uniprot_id, token, 
                    sizeof((*targets_ptr)[actual_count].uniprot_id) - 1);
            (*targets_ptr)[actual_count].uniprot_id[sizeof((*targets_ptr)[actual_count].uniprot_id) - 1] = '\0';
        }
        
        actual_count++;
    }
    
    fclose(file);
    printf("Successfully loaded %d target records from %s\n", actual_count, filename);
    return actual_count;
}


//Obtain target information based on composite name (automatically extract PDB ID)
struct target_info* get_target_info(const char *compound_name, struct target_info *targets, int target_count) {
    char extracted_pdb_id[20];
    extract_pdb_id(compound_name, extracted_pdb_id);
    
  //  printf("Looking for target: %s (extracted PDB ID: %s)\n", compound_name, extracted_pdb_id);
    
    int i;
    for (i = 0; i < target_count; i++) {
        if (strcmp(targets[i].pdb_id, extracted_pdb_id) == 0) {
            return &targets[i];
        }
    }
    return NULL;
}


 


 




















//Modify the function for printing detailed information of a single target 
void print_single_target_info(int rank, const char *compound_name, struct target_info *target) {
    char extracted_pdb_id[20];
    extract_pdb_id(compound_name, extracted_pdb_id);
    
    printf("\n%d. PDB ID: %s\n", rank, extracted_pdb_id);
    if (target != NULL) {
        printf("   Protein Names: %s\n", target->protein_names);
        printf("   Gene Names: %s\n", target->gene_names);
        //Add Uniprot ID display 
        if (target->uniprot_id[0] != '\0') {
            printf("   Uniprot ID: %s\n", target->uniprot_id);
        } else {
            printf("   Uniprot ID: Not available\n");
        }
    } else {
        printf("   Target Information: No information found\n");
    }
    printf("   %s\n", string(50, '-'));
}
//Print formatted separator lines 
char* string(int length, char character) {
    static char line[100];
    int i;
    for (i = 0; i < length && i < 99; i++) {
        line[i] = character;
    }
    line[i] = '\0';
    return line;
}


void free_targets_info(struct target_info *targets) {
    if (targets != NULL) {
        free(targets);
    }
}

//Add Uniprot ID column in CSV file save function 
void save_results_to_csv(struct docking_result *docking_results, int dock_num, 
                        struct target_info *targets, int target_count, 
                        const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not open file %s for writing\n", filename);
        return;
    }
    
    //Write CSV header - Add UniprotoID column 
    fprintf(file, "Rank,PDB_ID,Best_Score,Average_Score,Worst_Score,");
    fprintf(file, "Delta_Best,Delta_Avg,Delta_Worst,Composite_Delta,Similarity_Rank,");
    fprintf(file, "Protein_Names,Gene_Names,Uniprot_ID,Prediction\n");  // 添加Uniprot_ID
    
    int i;
    for (i = 0; i < dock_num; i++) {
        char extracted_pdb_id[20];
        extract_pdb_id(docking_results[i].receptor_name, extracted_pdb_id);
        struct target_info *target = get_target_info(docking_results[i].receptor_name, targets, target_count);
        
        //Write basic docking data 
        fprintf(file, "%d,%s,%.2f,%.2f,%.2f,", 
                i + 1, 
                extracted_pdb_id,
                docking_results[i].docking_scores.best_score,
                docking_results[i].docking_scores.average_score,
                docking_results[i].docking_scores.worst_score);
        
        fprintf(file, "%.2f,%.2f,%.2f,%.2f,%d,", 
                docking_results[i].delta_best,
                docking_results[i].delta_avg,
                docking_results[i].delta_worst,
                docking_results[i].delta_composite,
                docking_results[i].similarity_rank);
        
        //Write target information 
        if (target != NULL) {
            //Escaping commas in CSV 
            char protein_names_escaped[1000];
            char gene_names_escaped[500];
            char uniprot_id_escaped[100];
            
            strcpy(protein_names_escaped, target->protein_names);
            strcpy(gene_names_escaped, target->gene_names);
            strcpy(uniprot_id_escaped, target->uniprot_id);
            
            //Replace commas with semicolons to avoid CSV parsing issues 
            int k;
            for (k = 0; protein_names_escaped[k] != '\0'; k++) {
                if (protein_names_escaped[k] == ',') protein_names_escaped[k] = ';';
            }
            for (k = 0; gene_names_escaped[k] != '\0'; k++) {
                if (gene_names_escaped[k] == ',') gene_names_escaped[k] = ';';
            }
            for (k = 0; uniprot_id_escaped[k] != '\0'; k++) {
                if (uniprot_id_escaped[k] == ',') uniprot_id_escaped[k] = ';';
            }
            
            fprintf(file, "\"%s\",\"%s\",\"%s\",", protein_names_escaped, gene_names_escaped, uniprot_id_escaped);
        } else {
            fprintf(file, "\"No information found\",\"No information found\",\"No information found\",");
        }
        
        //Write prediction results
        if (docking_results[i].delta_best >= 999.0f) {
            fprintf(file, "DATA MISSING - cannot evaluate\n");
        } else if (docking_results[i].delta_best < -1.0) {
            fprintf(file, "HIGH potential target\n");
        } else if (docking_results[i].delta_best < -0.5) {
            fprintf(file, "GOOD potential target\n");
        } else if (docking_results[i].delta_best < 0.5) {
            fprintf(file, "MEDIUM potential target\n");
        } else if (docking_results[i].delta_best < 1.5) {
            fprintf(file, "LOW potential target\n");
        } else {
            fprintf(file, "VERY LOW potential target\n");
        }
    }
    
    fclose(file);
    printf("Results saved to %s\n", filename);
}
