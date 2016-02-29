#include "LoadConfig.h"

/* Function to read a file and creates a hash table.
The reading format is set to the model "key = value" line by line.*/
LoadConfig* file2map(const char* FileName)
{
	FILE *Archx = fopen(FileName, "r");

	LoadConfig *conf = Malloc(LoadConfig, 1);

	if(Archx == NULL)
    {
		fatal_error("Error reading file, function file2map. Check it please. \n");
        return NULL;
    }

	int lines = amount_line(Archx);

	conf->table = create_hash(lines);
	conf->getParameterChar = search_str_hash;

	char *key = Malloc(char, MAX_LINE_FILE);
	char *value = Malloc(char, MAX_LINE_FILE);

	for(int i = 0; i < conf->table->TABLE_SIZE; i++)
	{
		fscanf(Archx, "%s = %[^\n]s", key, value);
		insert_str_hash(conf->table, key, value);
	}

	fclose(Archx);
	free(key);
	free(value);

	return conf;
}

void close_conf(LoadConfig *conf)
{
	if(conf != NULL)
	{
		free_hash(conf->table);
		conf->getParameterChar = NULL;
		free(conf);
	}
	
}


