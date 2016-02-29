#ifndef MAPHASH_H
#define MAPHASH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "messages.h"

/*For performance reasons our hash table only stores the address of a generic type or struct. 
This avoids the excessive spending of memory and also allows you to choose the type of data to be stored .*/

struct s_HashItem
{
    void *hash_key;
    void *hash_item;
};
typedef struct s_HashItem HashItem_t;

struct s_HashTable
{
    int amount, TABLE_SIZE;
    HashItem_t **value;
};
typedef struct s_HashTable HashTable_t;

int key_string(const char *string);
int hashkey_divided(int key, int TABLE_SIZE);
int hashkey_multiplied(int key, int TABLE_SIZE);
int double_hash(int key, int TABLE_SIZE, int number_try, int Hash1);
HashTable_t *create_hash(int TABLE_SIZE);
void free_hash(HashTable_t *hash);
void insert_str_hash(HashTable_t *table, char *string_key, char *string_value);
char *search_str_hash(HashTable_t *table, const char *string);
int amount_line(FILE *arch);

#endif
