#include "maphash.h"

/*Calculating a key from a string , we have considered the case sensitive, otherwise anagrams could be considered as equal words*/
/*"value" variable is initialized with an integer, in the case "17" , only to the sum  doesn't always start with 0. Choose any integer, but give preference to prime numbers*/
/* The value 31 was chosen because it's an odd prime. Reference: hash function used in Java, for more details see the book: Effective Java by Joshua Bloch*/
/*Chapter 3, Item 9: Always override hashcode when you override equals, page 48*/

int key_string(const char *string)
{
	int value=17, str_size;

	str_size = strlen(string);

	for(int i = 0; i < str_size; i++)
	{
		value = 31 * value + (int) string[i];
	}

	return value;
}

/*Two hash functions for hashing. They're also used to treat collision ( Doublehash )*/
/*Modular function with Bitwise "AND" to remove the signal bit, preventing overflow"*/
/*Multiplicative hashing, aux is a real number between 0 and 1 (0 < aux < 1), in the case, Golden ratio - 1 =  0.6180339887*/

int hashkey_divided(int key, int TABLE_SIZE)
{
	return (key & 0xFFFFFFF) % TABLE_SIZE;
}

int hashkey_multiplied(int key, int TABLE_SIZE)
{
	float aux = 0.6180339887;
	float value = key * aux;
	value = value - (int) value;
	return (int) (value * TABLE_SIZE);
}

int double_hash(int key, int TABLE_SIZE, int number_try, int Hash1)
{
	int Hash2 = hashkey_divided(key, TABLE_SIZE-1) + 1;
	return ((Hash1 + number_try*Hash2) & 0xFFFFFFF) % TABLE_SIZE;
}

/* Create the hash table. 
Advice: I's very important to choose the size of the hash table. We leave it arbitrarily , but preferably the prime numbers in order that they reduce the likelihood of collisions. 
Avoid powers-of-two */

HashTable_t *create_hash(int TABLE_SIZE)
{
	HashTable_t *hash = Malloc(HashTable_t,1);

	if(hash != NULL)
	{
		hash->TABLE_SIZE = TABLE_SIZE;
		hash->amount = 0;
		hash->value = (HashItem_t**) malloc(TABLE_SIZE * sizeof(HashItem_t));

		if(hash->value == NULL)
		{
			free(hash);
			fatal_error("Allocation error in function create_hash/itens, please check it.\n");
			return NULL;
		}

		for(int i = 0; i < TABLE_SIZE; i++)
		{
			hash->value[i] = NULL;
		}

		return hash;
	}

	else
	{
		fatal_error("Allocation error in function create_hash, please check it.\n");
		return NULL;
	}
}

/*Destroys the hash table and all its elements*/

void free_hash(HashTable_t *hash)
{
	if(hash != NULL)
	{
		for(int i = 0; i < hash->TABLE_SIZE; i++)
		{
			if(hash->value[i] != NULL)
			{
				free(hash->value[i]->hash_key);
				free(hash->value[i]->hash_item);
				free(hash->value[i]);
			}
		}
		free(hash->value);
		free(hash);
	}
}

/* Insert and search function in the hash table with collision treatment. We're using strings as key and value, 
an integer or floating point could be used, simply modify the method of generating the key and the appropriate assignments of the elements*/

void insert_str_hash(HashTable_t *table, char *string_key, char *string_value)
{
	if(table == NULL || (table->amount == table->TABLE_SIZE))
	{
		fatal_error("Error in function insert_str_hash. Hash table doesn't exist or table's full, please check it.\n");
		return;
	}

	else
	{
		int position, newPosition;
		int key = key_string(string_key);

		position = hashkey_divided(key, table->TABLE_SIZE);

		for(int i = 0; i < table->TABLE_SIZE; i++)
		{
			newPosition = double_hash(position, table->TABLE_SIZE, i, hashkey_multiplied(key, table->TABLE_SIZE) );

			if(table->value[newPosition] == NULL)
			{
				HashItem_t *aux = Malloc(HashItem_t, 1);
				aux->hash_key = strdup(string_key);
				aux->hash_item = strdup(string_value);
				table->value[newPosition] = aux;
				table->amount++;
				return;

			}
		}

	}
}


char *search_str_hash(HashTable_t *table, const char *string)
{
	if(table == NULL)
	{
		fatal_error("Error in function search_str_hash. Hash table doesn't exist, please check it.\n");
		return NULL;
	}

	int position, newPosition;
	int key = key_string(string);

	position = hashkey_divided(key, table->TABLE_SIZE);

	for(int i = 0; i < table->TABLE_SIZE; i++)
	{
		newPosition = double_hash(position, table->TABLE_SIZE, i, hashkey_multiplied(key, table->TABLE_SIZE) );

		if(strcmp((char*)table->value[newPosition]->hash_key, string) == 0)
		{
			return (char*)table->value[newPosition]->hash_item;
		}

	}

	return (char*) "\0";
}

/*Function to count lines of a file. We use it to set the size of our hash table, but as stated above, the choice is arbitrary.*/

int amount_line(FILE *arch)
{
	int count = 0;
	char c;

	while((c = getc(arch)) != EOF)
	{
		if(c == '\n')
			count++;
	}
	rewind(arch);
	return count+1;
}
