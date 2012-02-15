
%module (docstring="A wrapper for Venice") venicepy
%{
#include "main.h"
%}

%include "main.h"

extern int main(int argc, char **argv);
extern void testPython();

