#include<iostream>
#include<iomanip>

#include<fftw3.h>

using namespace std;
int main(){
    int N=10;

    fftw_complex *in,*out;
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    if((in==NULL)||(out==NULL)){
        printf("Error:insufficient available memory\n");
    }
    else{
        for(int i=0; i<N; i++){
            in[i][0] = i+1;
            in[i][1] = 0;
        }
    }

    fftw_plan  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_execute(p); /* repeat as needed */
    fftw_destroy_plan(p);
    fftw_cleanup();

    for(int i=0;i<N;i++){/*OUTPUT*/
        cout<<setprecision(6)<<setiosflags(ios::fixed);
        cout<<in[i][0]<<", "<<in[i][1]<<"i"<<endl;
    }
    cout<<endl;
    for(int i=0;i<N;i++){/*OUTPUT*/
        cout<<setprecision(6)<<setiosflags(ios::fixed);
        cout<<out[i][0]<<", "<<out[i][1]<<"i"<<endl;
    }

    if(in!=NULL) fftw_free(in);
    if(out!=NULL) fftw_free(out);

    return 0;
}
