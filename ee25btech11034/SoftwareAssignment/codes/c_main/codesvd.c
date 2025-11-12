#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double* transpose(int m,int n,double*arr){
    double*arr_T=(double*)malloc(n*m*sizeof(double));
    if(arr_T==NULL){
        printf("Error: Failed to allocate memory for U sorting.\n");
        exit(1);
    }
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
           arr_T[j*m+i]=arr[i*n+j];
        }
    }
    return arr_T;
}
//old multiply which is slower
// double* multiply(int m,int k,int n,double*arr1,double*arr2){
///     double*arr3=(double*)malloc(m*n*sizeof(double));
//     if(arr3==NULL){
//         printf("Error: Failed to allocate memory for U sorting.\n");
//         exit(1);
//     }
//       for(int i=0;i<m;i++){
//         for(int j=0;j<n;j++){
//            arr3[i*n+j]=0;
//            for(int a=0;a<k;a++){
//                arr3[i*n+j]+=arr1[i*k+a]*arr2[a*n+j];
//            }
//         }
//       }
//       return arr3;
// }
double* multiply(int m,int k,int n,double*arr1,double*arr2){
    double*arr3=(double*)malloc(m*n*sizeof(double));
    if(arr3==NULL){
        printf("Error: Failed to allocate memory for arr3.\n");
        exit(1);
    }
    int tile=64;


    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            arr3[i*n+j]=0;
        }
    }
    for(int i1=0;i1<m;i1+=tile){
        for(int j1=0;j1<n;j1+=tile){
            for(int a1=0;a1<k;a1+=tile){
                for(int i=i1;i<i1+tile&&i<m;i++){
                    for(int j=j1;j<j1+tile&&j<n;j++){
                        for(int a=a1;a<a1+tile&&a<k;a++){
                            arr3[i*n+j]+=arr1[i*k+a]*arr2[a*n+j];
                        }
                    }
                }
            }
        }
    }
    return arr3;
}

void merge(double*values,int*indices,int left,int mid,int right){
    int i,j,k;
    int n1=mid-left+1;
    int n2=right-mid;

    // Create temporary arrays for values and indices
    double*l_val=(double*)malloc(n1*sizeof(double));
    int*l_idx=(int*)malloc(n1*sizeof(int));
    double*r_val=(double*)malloc(n2*sizeof(double));
    int*r_idx=(int*)malloc(n2*sizeof(int));

    if(l_val==NULL||l_idx==NULL||r_val==NULL||r_idx==NULL){
        printf("Error: Failed to allocate memory in merge.\n");
        exit(1);
    }

    // Copy data to temp arrays
    for(i=0;i<n1;i++){
        l_val[i]=values[left+i];
        l_idx[i]=indices[left+i];
    }
    for(j=0;j<n2;j++){
        r_val[j]=values[mid+1+j];
        r_idx[j]=indices[mid+1+j];
    }

    // Merge the temp arrays back into values[] and indices[]
    i=0; // Initial index of first subarray
    j=0; // Initial index of second subarray
    k=left; // Initial index of merged subarray
    while(i<n1&&j<n2){
        // Sort in descending order
        if(l_val[i]>=r_val[j]){
            values[k]=l_val[i];
            indices[k]=l_idx[i];
            i++;
        }else{
            values[k]=r_val[j];
            indices[k]=r_idx[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of L[]
    while(i<n1){
        values[k]=l_val[i];
        indices[k]=l_idx[i];
        i++;
        k++;
    }

    // Copy the remaining elements of R[]
    while(j<n2){
        values[k]=r_val[j];
        indices[k]=r_idx[j];
        j++;
        k++;
    }
    free(l_val);
    free(l_idx);
    free(r_val);
    free(r_idx);
}

// Sorts an array of values and a parallel array of indices using Merge Sort.
void mergeSort(double*values,int*indices,int left,int right){
    if(left<right){
        int mid=left+(right-left)/2;

        // Sort first and second halves
        mergeSort(values,indices,left,mid);
        mergeSort(values,indices,mid+1,right);

        merge(values,indices,left,mid,right);
    }
}
void print_matrix(int m,int n,double*arr,const char*label){
    printf("%s (%d x %d):\n",label,m,n);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            printf("%10.6f ",arr[i*n+j]);
        }
        printf("\n");
    }
    printf("\n");
}
int min(int a,int b)
{
    return(a<b)?a:b;
}
void jacobi(int n,double*arr,double*V)
{
    double theta;
    double tan2;
    double t;
    double c;
    double s;
    double EPSILON=1.0e-10;
    double old[n];
    for(int i=0;i<n;i++){
        old[i]=arr[n*i+i];
    }

for(int steps=0;steps<100;steps++)
    {
        for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            if(fabs(arr[n*i+j])<EPSILON){
                continue;
            }
            if(fabs(arr[n*i+i]-arr[n*j+j])<EPSILON)
            {
                t=1.0;
                theta=(arr[n*i+j]>0)?acos(-1)/4:-acos(-1)/4;
            }
            else{
                theta=0.5*atan2(2*arr[n*i+j],(arr[n*i+i]-arr[n*j+j]));
                t=tan(theta);
            }
            c=cos(theta);
            s=sin(theta);
            double arr_ii=arr[n*i+i];
            double arr_jj=arr[n*j+j];
            double arr_ij=arr[n*i+j];
            arr[n*i+i]=c*c*arr_ii+2.0*s*c*arr_ij+s*s*arr_jj;
            arr[n*j+j]=s*s*arr_ii-2.0*s*c*arr_ij+c*c*arr_jj;
            //arr[i][j] = sin*cos*(arr_ii-arr_jj)+(cos*cos-sin*sin)*arr_ij;
            arr[n*i+j]=0;
            arr[n*j+i]=arr[n*i+j];
            for(int k=0;k<n;k++)
            {
                if(k!=i&&k!=j){
                    double arr_ki=arr[k*n+i];
                    double arr_kj=arr[k*n+j];
                    arr[k*n+i]=c*arr_ki+s*arr_kj;
                    arr[k*n+j]=s*arr_ki-c*arr_kj;
                    arr[i*n+k]=arr[k*n+i];
                    arr[j*n+k]=arr[k*n+j];
                }
            }
            for(int k=0;k<n;k++){
                double v_ki=V[k*n+i];
                double v_kj=V[k*n+j];
                V[k*n+i]=c*v_ki+s*v_kj;
                V[k*n+j]=-s*v_ki+c*v_kj;
            }

        }
    }
    double disum=0.0;
    for(int i=0;i<n;i++){
        disum+=fabs(arr[n*i+i]-old[i]);
        old[i]=arr[n*i+i];
    }
    if(disum<EPSILON){
        break;
    }
    }
    // for(int i=0;i<n;i++){
    //         for(int j=0;j<n;j++){
    //                printf("%lf ",arr[n*i+j]);
    //         }
    //         printf("\n");
    // }
}
void svd(int m,int n,double*A,double*U,double*S,double*V){
    // compute A_T (n x m)
    double*A_T=transpose(m,n,A);

    // compute A_T * A (n x n, symmetric).
    double*ATA=multiply(n,m,n,A_T,A);

    // V starts as identity (caller allocated V)
    for(int i=0;i<n;i++){
        V[i*n+i]=1.0;
    }

    // calculate V and Eigenvalues from A_T * A
    jacobi(n,ATA,V);

    // sort V based on eigenvalues (diagonal of ATA)
    double*v_values=(double*)malloc(n*sizeof(double));
    int*v_indices=(int*)malloc(n*sizeof(int));
    double*V_temp=(double*)malloc(n*n*sizeof(double));
    if(v_values==NULL||v_indices==NULL||V_temp==NULL){
        printf("Error: Failed to allocate memory for V sorting.\n");
        exit(1);
    }
    for(int i=0;i<n*n;i++){
        V_temp[i]=V[i];
    }

    for(int i=0;i<n;i++){
        v_values[i]=ATA[i*n+i];
        v_indices[i]=i;
    }
    mergeSort(v_values,v_indices,0,n-1);

    // Reorder V columns
    for(int i=0;i<n;i++){ // i is the new column index (sorted)
        int old_col=v_indices[i];
        for(int row=0;row<n;row++){
            V[row*n+i]=V_temp[row*n+old_col];
        }
    }
    free(V_temp);
    
    //  v contains the eigenvectors of A_T*A (which is V), sorted by eigenvalue
    //  v_values contains the sorted eigenvalues (sigma^2)

    // Compute A * A_T (m x m, symmetric)
    // We need a copy because jacobi modifies the input.
    double*AAT=multiply(m,n,m,A,A_T);

    // U starts as identity 
    for(int i=0;i<m;i++){
        U[i*m+i]=1.0;
    }

    // calculate U 
    //jacobi(m, AAT, U);
    double*AV=multiply(m,n,n,A,V);

    // calculate U columns: u_i = (1 / s_i) * (A*V)_i
    for(int i=0;i<min(m,n);i++){
        // use the sorted eigenvalues from v_values
        double singular_value=sqrt(fabs(v_values[i]));
        S[i*n+i]=singular_value;
    }
    for(int j=0;j<min(m,n);j++){ 
        double s_j=S[j*n+j];
        double s_inv=1.0;

        if(fabs(s_j)>1.0e-12){ 
            s_inv=1.0/s_j;
        }

        // copy the j-th column of AV into the j-th column of U scale by 1/s_j
        for(int i=0;i<m;i++){ 
            U[i*m+j]=AV[i*n+j]*s_inv;
        }
    }

    // U contains the eigenvectors of A*A_T (which is U), sorted by eigenvalue

    // ceate S (m x n) from the sorted eigenvalues (diagonal of ATA)
    int k=(m<n)?m:n;

    // correct sign
    double*US=multiply(m,m,n,U,S);

    for(int i=0;i<min(m,n);i++){ // For each column
        double dot_product=0.0;
        for(int row=0;row<m;row++){
            // Dot product of column i of AV and column i of US
            dot_product+=AV[row*n+i]*US[row*n+i];
        }

        if(dot_product<0.0){
            // Flip the sign of the i-th column of U
            for(int row=0;row<m;row++){
                U[row*m+i]*=-1.0;
            }
        }
    }

    free(A_T);
    free(ATA);
    free(AAT);
    free(AV);
    free(US);
    free(v_values);
    free(v_indices);

}

double* topk(int m,int n,int k,double*U,double*V,double*S,double*U1,double*V1T,double*S1){
    for(int i=0;i<m;i++){
        for(int j=0;j<k;j++){
            U1[i*k+j]=U[i*m+j];
        }
    }
    for(int i=0;i<k;i++){
        for(int j=0;j<n;j++){
            V1T[i*n+j]=V[j*n+i];
        }
    }
    for(int i=0;i<k;i++){
            S1[i*k+i]=S[i*n+i];
    }
    double*A1=multiply(m,k,k,U1,S1);
    double*A2=multiply(m,k,n,A1,V1T);
    free(A1);
    return A2;
}
double* read_image(const char*filename,int*width,int*height){
    int channels;
    // force grayscale by requesting 1 channel
    unsigned char*img_data=stbi_load(filename,width,height,&channels,1);

    if(img_data==NULL){
        printf("Error: Failed to load image '%s'\n",filename);
        printf("Reason: %s\n",stbi_failure_reason());
        return NULL;
    }
    printf("Image loaded successfully: %s\n",filename);
    printf("Dimensions: %d x %d, Original channels: %d\n",*width,*height,channels);

    int tot_pixels=(*width)*(*height);
    double*pixel_matrix=(double*)malloc(tot_pixels*sizeof(double));
    if(pixel_matrix==NULL){
        printf("Error: Failed to allocate memory for pixel matrix\n");
        stbi_image_free(img_data);
        return NULL;
    }

    // convert unsigned char (0-255) to double (0.0-255.0)
    // keep the range 0-255 to maintain precision
    for(int i=0;i<tot_pixels;i++){
        pixel_matrix[i]=(double)img_data[i];
    }

    stbi_image_free(img_data);
    printf("Image converted to double matrix successfully\n");
    return pixel_matrix;
}
int write_image(const char*filename,int width,int height,double*pixel_matrix){
    int tot_pixels=width*height;

    // allocate unsigned char array
    unsigned char*img_data=(unsigned char*)malloc(tot_pixels*sizeof(unsigned char));
    if(img_data==NULL){
        printf("Error: Failed to allocate memory for image writing\n");
        return 0;
    }

    // convert double to unsigned char with clamping
    for(int i=0;i<tot_pixels;i++){
        double pixel_val=pixel_matrix[i];

        // Clamp to [0, 255] range
        if(pixel_val<0.0){
            pixel_val=0.0;
        }else if(pixel_val>255.0){
            pixel_val=255.0;
        }

        // round and convert to unsigned char
        img_data[i]=(unsigned char)(pixel_val+0.5);
    }

    // get file format from extension
    int result=0;
    const char*ext=strrchr(filename,'.');

    if(ext==NULL){
        printf("Error: No file extension provided. Using PNG format.\n");
        result=stbi_write_png(filename,width,height,1,img_data,width);
    }else if(strcmp(ext,".png")==0||strcmp(ext,".PNG")==0){
        result=stbi_write_png(filename,width,height,1,img_data,width);
    }
    else if(strcmp(ext,".jpg")==0||strcmp(ext,".JPG")==0||strcmp(ext,".jpeg")==0||strcmp(ext,".JPEG")==0){
        result=stbi_write_jpg(filename,width,height,1,img_data,90); // 90 = quality
    }
    else{
        printf("Warning: Unknown extension '%s'. Using PNG format.\n",ext);
        result=stbi_write_png(filename,width,height,1,img_data,width);
    }

    free(img_data);

    if(result==0){
        printf("Error: Failed to write image '%s'\n",filename);
        return 0;
    }

    printf("Image saved successfully: %s\n",filename);
    return 1;
}

int main(){
    char ip_fname[1024];
    char op_fname[1024];
    int k;
    printf("Enter the input image path: ");
    scanf("%s",ip_fname);

    printf("Enter output image path: ");
    fflush(stdout); // Force the prompt to display
    scanf("%1023s", op_fname);

    printf("Enter k value: ");
    fflush(stdout); // Force the prompt to display
    scanf("%d",&k);


    if(k<=0){
        printf("Error: k value must be positive\n");
        return 1;
    }

    // read the image
    int width,height;
    double*pixel_matrix=read_image(ip_fname,&width,&height);

    if(pixel_matrix==NULL){
        return 1;
    }

    printf("\n=== Starting SVD Compression ===\n");
    printf("Original dimensions: %d x %d\n",height,width);
    printf("Compression rank: k = %d\n",k);

    // check if k is valid
    int max_k=(height<width)?height:width;
    if(k>max_k){
        printf("Warning: k = %d is larger than min(m,n) = %d. Setting k = %d\n",k,max_k,max_k);
        k=max_k;
    }

    // allocate SVD matrices
    int m=height;
    int n=width;

    double*U=(double*)calloc(m*m,sizeof(double));
    double*S=(double*)calloc(m*n,sizeof(double));
    double*V=(double*)calloc(n*n,sizeof(double));
    double*U1=(double*)calloc(m*k,sizeof(double));
    double*S1=(double*)calloc(k*k,sizeof(double));
    double*V1T=(double*)calloc(k*n,sizeof(double));

    if(U==NULL||S==NULL||V==NULL||U1==NULL||S1==NULL||V1T==NULL){
        printf("Error: Failed to allocate memory for SVD matrices\n");
        free(pixel_matrix);
        
        free(U);
        free(S);
        free(V);
        free(U1);
        free(S1);
        free(V1T);
        return 1;
    }

    // perform SVD
    printf("Performing SVD decomposition...\n");
    svd(m,n,pixel_matrix,U,S,V);

    // reconstruct with top-k singular values
    printf("Reconstructing image with top-%d singular values...\n",k);
    double*compressed_matrix=topk(m,n,k,U,V,S,U1,V1T,S1);

    // write the compressed image
    printf("\n=== Writing Output Image ===\n");
    int success=write_image(op_fname,width,height,compressed_matrix);
    double tot_energy_sq=0.0;
    double err_energy_sq=0.0;
    for(int i=0;i<max_k;i++){
        double s_val=S[i*n+i];
        double s_val_sq=s_val*s_val;
        tot_energy_sq+=s_val_sq;
        if(i>=k){
            err_energy_sq+=s_val_sq;
        }
    }
    double frob_norm_error=sqrt(err_energy_sq);
    double frob_norm_total=sqrt(tot_energy_sq);
    double percent_retained=0.0;
    if(frob_norm_total>0){
            percent_retained=100.0*(1.0-(err_energy_sq/tot_energy_sq));
    }

    if(success){
        printf("Input: %s\n",ip_fname);
        printf("Output: %s\n",op_fname);
        printf("Rank: %d / %d (%.2f%% of singular values retained)\n",k,max_k,(100.0*k)/max_k);
        printf("Frobenius Norm of Total Image: %12.2f\n",frob_norm_total);
        printf("Frobenius Norm of Error (A-A2): %12.2f\n",frob_norm_error);

    }

    free(pixel_matrix);
    free(U);
    free(S);
    free(V);
    free(U1);
    free(S1);
    free(V1T);
    free(compressed_matrix);

    return success?0:1;
}