#ifndef KAHN_INCLUDED
#define KAHN_INCLUDED

/*      
        Interface implémentation minimale permettant de construire des réseaux
        de processus de Kahn.

        Les processus sont représenté par la structure process comportant 2 champs:
            - la fonction associée au processus
            - les paramètres à donner à la fonction lors de son éxécution
            
        Les objets mis à disposition par l'interface sont:
            - Kahn_Datatype // types mis à disposition pour les communications
            - channel       // canaux de communication entre processus (pipes)
            - new_channel() // fonction initiation d'un canal
            - put()         // envoi sur un canal
            - get()         // réception depuis un canal
            - doco()        // éxécution parallèle d'une liste de processus
            - run()         // éxécute un processus
            - return()      // renvoi un processus renvoyant une donnée et terminant
            - bind()        // bind 2 processus et renvoi le processus résultant
*/

typedef enum { 
    KAHN_INT, 
    KAHN_FLOAT,
    KAHN_DOUBLE,
    _nb_kahn_datatype,
} Kahn_Datatype;

typedef void *(*process_function)();

typedef struct {
    process_function f;
    void *arg;
} process;

typedef struct {
    int fd_in;
    int fd_out;
} channel;

channel *new_channel();

void put(void *value, int cnt, channel *chan_out, Kahn_Datatype dtype);
void get(void *value, int cnt, channel *chan_in, Kahn_Datatype dtype);

void doco(int nb_proc, process *processes[]);

process *return_(void *value);
process *bind(process *p1, process *p2);
void *run(process *p);

#endif
