#include <stdio.h>
#include <time.h>
#include <mpi.h>
#define WIDTH 640
#define HEIGHT 480
#define MAX_ITER 255

struct complex{
  double real;
  double imag;
};
int cal_pixel(struct complex c) {
    double z_real = 0;
    double z_imag = 0;

    double z_real2, z_imag2, lengthsq;

    int iter = 0;
    do {
        z_real2 = z_real * z_real;
        z_imag2 = z_imag * z_imag;

        z_imag = 2 * z_real * z_imag + c.imag;
        z_real = z_real2 - z_imag2 + c.real;
        lengthsq =  z_real2 + z_imag2;
        iter++;
    }
    while ((iter < MAX_ITER) && (lengthsq < 4.0));

        return iter;

}
void save_pgm(const char *filename, int image[HEIGHT][WIDTH]) {
    FILE* pgmimg; 
    int temp;
    pgmimg = fopen(filename, "wb"); 
    fprintf(pgmimg, "P2\n"); // Writing Magic Number to the File   
    fprintf(pgmimg, "%d %d\n", WIDTH, HEIGHT);  // Writing Width and Height
    fprintf(pgmimg, "255\n");  // Writing the maximum gray value 
    int count = 0; 
    
    for (int i = 0; i < HEIGHT; i++) { 
        for (int j = 0; j < WIDTH; j++) { 
            temp = image[i][j]; 
            fprintf(pgmimg, "%d ", temp); // Writing the gray values in the 2D array to the file 
        } 
        fprintf(pgmimg, "\n"); 
    } 
    fclose(pgmimg); 
} 

int main(int argc, char const *argv[])
{
    int image[HEIGHT][WIDTH];
    struct complex c;
	MPI_Init(NULL,NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int i,j,l,r,row_recieved,s=HEIGHT/(world_size-1),slave;
    int row[WIDTH+2];
    int data_tag=0, result_tag=0, source_tag=0,terminator_tag=-1;
    double time;
    if(world_rank==0){
        clock_t start_time = clock();
        int count=0;
        r=0;
        for(i=1;i<world_size;i++){
            MPI_Send(&r,sizeof(int),MPI_INT,i,data_tag,MPI_COMM_WORLD);
            count++;
            r++;
        }
        do{
            MPI_Recv(&row,(WIDTH+2)*sizeof(int),MPI_INT,MPI_ANY_SOURCE,result_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            slave=row[WIDTH+1];
            row_recieved=row[WIDTH+2];
            count--;
            if(r<HEIGHT){
                MPI_Send(&r,sizeof(int),MPI_INT,slave,data_tag,MPI_COMM_WORLD);
                r++;
                count++;
            }
            else{
                MPI_Send(&r,sizeof(int),MPI_INT,slave,terminator_tag,MPI_COMM_WORLD);
            }
            for(int k=0;k<WIDTH;k++){
                image[row_recieved][k]=row[k];
            }
        }while(count>0);
        clock_t end_time = clock();
        save_pgm("mandelbrot_static.pgm", image);
        time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
        printf("The average execution time is: %f ms", time*1000);
    }
    else{
        MPI_Recv(&r,sizeof(int),MPI_INT,0,source_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        while(r>=0){
            c.imag=(r-HEIGHT/2.0)*4.0/HEIGHT;
            for(i=0;i<WIDTH;i++){
                c.real=(i-HEIGHT/2.0)*4.0/WIDTH;
                row[i]=cal_pixel(c);
            }
            row[WIDTH]=world_rank;
            row[WIDTH+1]=r;
            MPI_Send(&row,(WIDTH+2)*sizeof(int),MPI_INT,0,result_tag,MPI_COMM_WORLD);
            MPI_Recv(&r,sizeof(int),MPI_INT,0,source_tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    MPI_Finalize();

	return 0;
}