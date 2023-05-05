#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <string.h>
#include <assert.h>

#include "kahn.h"

/*      Implémentation de l'interface kahn.h
        Objets encapsulés:
            - _kahn_data_size   // donne les tailles en octets des types utilisés
            - _pid_list         // liste des pid des processus forkés par doco
            - _return_process() // fonction du processus retourné par return_()

        DEBUG: variable de préprocesseur contrôlant le niveau de verbose
*/


//#define DEBUG

/* ***   COMMUNICATION   *** */

channel *new_channel() {
    channel *chan = (channel*) malloc(sizeof(channel));
    int fd[2];
    pipe(fd);   // initiation d'une pipe, retourne les file descriptor de in et out dans fd
    chan->fd_in = fd[0];
    chan->fd_out = fd[1];
#ifdef DEBUG
    printf("KAHN\t[%d] NEW CHANNEL ; fd_in = %d\tfd_out = %d\n", (int)getpid(), fd[1], fd[0]);
#endif
    return chan;
    }

static int _kahn_data_size[] = {
    sizeof(int),
    sizeof(float),
    sizeof(double),
    };
    
void put(void *value, int cnt, channel *chan, Kahn_Datatype dtype) 
{
#ifdef DEBUG
    printf("KAHN\t[%d] Entering PUT ; fd_out = %d; cnt = %d\n", (int)getpid(), chan->fd_out, cnt);

#endif
    assert((0 <= dtype) && (dtype < _nb_kahn_datatype));

    int total_write = 0;
    total_write = write(chan->fd_out, value,cnt* _kahn_data_size[dtype]);
    total_write = total_write;
#ifdef DEBUG
    printf("KAHN\t[%d] Exiting PUT ; cnt_write = %d\n", (int)getpid(), total_write/_kahn_data_size[dtype]);
#endif
}

void get(void *value, int cnt,channel *chan, Kahn_Datatype dtype) 
/*
    Buffered read
*/
{
#ifdef DEBUG
    printf("KAHN\t[%d] Entering GET ; fd_in = %d ; cnt = %d\n", (int)getpid(),chan->fd_in, cnt);
#endif
    assert((0 <= dtype) && (dtype < _nb_kahn_datatype));

    int total_read = 0;

    while(total_read < cnt* _kahn_data_size[dtype]) {
        total_read += read(chan->fd_in, value+(total_read/_kahn_data_size[dtype]), cnt* _kahn_data_size[dtype] - total_read);
    }


#ifdef DEBUG
    printf("KAHN\t[%d] Exiting GET ; cnt_read = %d\n", (int)getpid(), total_read/_kahn_data_size[dtype]);
#endif
}

/* ***   RETURN BIND RUN *** */

void *_return_process(void *value)
{
    return value;
}

process *return_(void *value)
{
    process *p = malloc(sizeof(*p));
    p->f = _return_process;
    p->arg = value;
    return p;
}

void *run(process *p)
{
    return (p->f)(p->arg);
}

process *bind(process *p1, process *p2)
{
    process *p = malloc(sizeof(*p));
    p->f = p2->f;
    p->arg = run(p1);
    return p;
}


/* ***   DOCO   *** */

int _kahn_is_main = 1;

void doco(int nb_proc, process *processes[]) {

#ifdef DEBUG
    printf("KAHN\t[%d] Entering doco, nb_proc = %d\n", (int)getpid(),nb_proc);
    fflush(stdout);
#endif

    if (nb_proc > 0) {
        pid_t pid = fork();
        if (pid == 0) {
            /* SON */
            _kahn_is_main = !_kahn_is_main;
            run(*processes);
        }
        else {
            /* FATHER */
#ifdef DEBUG
            printf("KAHN\t[%d] Creating [%d]\n", (int)getpid(), (int)pid);
#endif
            doco(--nb_proc, ++processes);

            /* Seul le processus initial attend ses enfants */
            if (!_kahn_is_main) {return; }

            wait(NULL);
#ifdef DEBUG
            printf("KAHN\t[%d] Something has ended\n", (int)getpid());
#endif
        }
    }
}