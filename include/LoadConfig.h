#ifndef LOADCONFIG_H
#define LOADCONFIG_H

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "maphash.h"
#include "messages.h"

struct s_conf
{	
	char* (*getParameterChar) (HashTable_t *table, const char *string);
	HashTable_t *table;
};
typedef struct s_conf LoadConfig;

LoadConfig* file2map(const char* FileName);
void close_conf(LoadConfig *conf);

#endif