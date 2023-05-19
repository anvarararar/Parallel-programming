#include <stdio.h>
#include <malloc.h>
#include <mpi.h>

extern double func  (double t, double x);
extern double fi    (double x);
extern double ksi   (double t);

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;

typedef struct 
{
    double x_start;
    size_t x_size, t_size;
    double* array;
} TaskInfo;

TaskInfo* TaskCreator (int rank, int size)
{
    TaskInfo* task = (TaskInfo*) calloc(1, sizeof (*task));
    if (!task)
        return NULL;

    task->x_start = x_max * rank / size; // начало участка нашего потока
    // размеры сетки:
    size_t x_size = x_max / x_step;
    task->t_size = t_max / t_step;
    
    task->x_size = x_size / size; // только наш участок
    if (rank == size - 1)
        task->x_size += x_size % size;

    task->array = (double*) calloc (task->x_size * task->t_size, sizeof (double));
    if (!task->array)
    {
        free (task);
        return NULL;
    }
    return task;
}

void ArrayInitialize (TaskInfo* task, int rank)
{
    double x_now = task->x_start;
    for (size_t i_x = 0; i_x < task->x_size; i_x++, x_now += x_step)
        task->array[i_x] = fi(x_now);
    if (rank == 0)
    {
        double t_now = 0.0;
        for (size_t i_t = 0; i_t < task->t_size; i_t++, t_now += t_step)
            task->array[i_t * task->x_size] = ksi(t_now);
    }
}

void TaskDestroyer (TaskInfo* task)
{
    if (task)
    {
        free (task->array);
        free (task);
    }
}

void CountByLeftAngle (int rank, int size, TaskInfo* task)
{
    double t_now = t_step;
    int first = !rank,
        not_last = (rank < size - 1);

    for (size_t t_i = 1; t_i < task->t_size; t_i++, t_now += t_step)
    {
        double u_prev = task->array[(t_i - 1) * task->x_size];
        if (!first)
        {
            MPI_Recv (&u_prev, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (not_last)
            MPI_Send (&(task->array[task->x_size * t_i - 1]), 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

        double x_now = task->x_start + x_step * first;
        for (size_t x_i = first; x_i < task->x_size; x_i++, x_now += x_step)
        {
            double u_now = task->array[x_i + task->x_size * (t_i - 1)];
            task->array[x_i + task->x_size * t_i] = 
                u_now + t_step * (func (t_now, x_now) - (u_now - u_prev) / x_step);
            u_prev = u_now;
        }
    }
}

void PrintToFile (TaskInfo* task, int rank, int size)
{
    // Записываем в csv последовательно из всех потоков:
    int wait = 0;
    if (rank != 0)
    {
        MPI_Recv (&wait, 1, MPI_INT, rank - 1, MPI_ANY_TAG,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    const char* flag = rank ? "a" : "w";
    FILE *fp = fopen("output.csv", flag);
    if (!fp)
    {
        printf ("Unable to open file.\n");
        return;
    }

    if (rank == 0)
        fprintf(fp, "x,t,u\n"); // записываем заголовок

    size_t t_print_step = task->t_size > 100 ? task->t_size / 100 : 1,
           x_print_step = task->x_size > 100 ? task->x_size / 100 : 1;

    for (size_t t_i = 0; t_i < task->t_size; t_i += t_print_step)
    {
        for (size_t x_i = 0; x_i < task->x_size; x_i += x_print_step)
            fprintf(fp, "%.6f,%.6f,%.6f\n", task->x_start + x_i * x_step, t_i * t_step,
                    task->array[x_i + task->x_size * t_i]);
    }
    fclose(fp);

    if (rank < size - 1)
        MPI_Send (&wait, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
}

int main( int argc, char **argv )
{
    MPI_Init (&argc, &argv);
    int size = 0, rank = 0;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    TaskInfo* task = TaskCreator (rank, size);
    if (!task)
    {
        printf ("Unable to allocate memory.\n");
        MPI_Finalize();
        return 0;
    }

    double t_start, t_end;
    if(rank == 0)
    	t_start = MPI_Wtime();

    ArrayInitialize (task, rank);
    CountByLeftAngle (rank, size, task);
    
    if(rank == 0)
    {
    	t_end = MPI_Wtime();
        printf ("Working time: %lf\n", t_end - t_start);
    }
    PrintToFile (task, rank, size);

    TaskDestroyer (task);
    MPI_Finalize();
    return 0;
}
