#include <stdlib.h>
#include <stdio.h>

typedef struct node{
    long L;
} node;

void change(node* tree){
    for(int i=0;i<10;i++){
        node curr=tree[i];
        curr.L=1;
    }
}

int main(){
    node* tree=(node*)malloc(10*sizeof(node));
    for(int i=0;i<10;i++)
        tree[i].L=0;
    change(tree);
    for(int i=0;i<10;i++)
        printf("%ld\n",tree[i].L);
}