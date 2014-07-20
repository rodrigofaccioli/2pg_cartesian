#include<stdio.h>
#include<stdlib.h>
#include <signal.h>

#include"messages.h"

_2PG_CARTESIAN_EXPORT
void fatal_error(const char *__restrict __message){
	char msg[10240];
	sprintf(msg,"\n%s \n",__message);
	fprintf(stderr,"%s\n",msg);
	raise(SIGABRT);
}

_2PG_CARTESIAN_EXPORT
void display_msg(const char *__restrict __message){
	char msg[10240];
	sprintf(msg,"%s",__message);
	fprintf(stdout,"%s",msg);
}
