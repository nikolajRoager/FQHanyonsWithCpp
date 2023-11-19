#include "mpo.hpp"
#include "approx.hpp"
#include <thread>

//n operator
mpo::mpo(
    uint64_t Lx,
    uint64_t Ly)
{
    is_n=true;

    L=Lx*Ly;
    d=2;
    dTwo=d*d;
    Dw=2;

    matrices=vector<sp_cx_mat>(L*dTwo,sp_cx_mat(Dw,Dw));//There are L*(d**2) matrices

    column_lookup=vector<uint>(1,0);//Not really worth anything here

    matrices[0]=sp_cx_mat(1,Dw);
    matrices[1]=sp_cx_mat(1,Dw);
    matrices[2]=sp_cx_mat(1,Dw);
    matrices[3]=sp_cx_mat(1,Dw);

    matrices[L*dTwo-4]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-3]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-2]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-1]=sp_cx_mat(Dw,1);


    //sig=0, sig1=0
    matrices[0](0,0)={0.0,0.0};
    matrices[0](0,1)={1.0,0.0};

    //Will be modified when we take the expectation value
    matrices[3](0,0) = {0.0,0.0};//sig=1, sig1=1
    matrices[3](0,1) = {1.0,0.0};//sig=1, sig1=1

    for (uint j = 1; j+1 < L; ++j)
    {

        matrices[0+dTwo*j](0,0)={1.0,0.0};
        matrices[0+dTwo*j](1,0)={0.0,0.0};matrices[0+dTwo*j](1,1)={1.0,0.0};

        matrices[3+dTwo*j](0,0)={1.0,0.0};
    //Will be modified when we take the expectation value
        matrices[3+dTwo*j](1,0)={0.0,0.0};matrices[3+dTwo*j](1,1)={1.0,0.0};
    }

    matrices[L*dTwo-4](0,0)={1.0,0.0};
    matrices[L*dTwo-4](1,0)={0.0,0.0};
    matrices[L*dTwo-1](0,0)={1.0,0.0};
    matrices[L*dTwo-1](1,0)={0.0,0.0};
}



//1D hopping
mpo::mpo(
    uint64_t _L,
    const vector<double>& V,
    complex<double> t)
{
    L=_L;
    d=2;
    dTwo=d*d;
    Dw=4;

    matrices=vector<sp_cx_mat>(L*dTwo,sp_cx_mat(Dw,Dw));//There are L*(d**2) matrices

    matrices[0]=sp_cx_mat(1,Dw);
    matrices[1]=sp_cx_mat(1,Dw);
    matrices[2]=sp_cx_mat(1,Dw);
    matrices[3]=sp_cx_mat(1,Dw);

    matrices[L*dTwo-4]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-3]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-2]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-1]=sp_cx_mat(Dw,1);


    column_lookup=vector<uint>(Dw-1,0);//Not really worth anything here



    //sig=0, sig1=0
    matrices[0](0,0)={V[0]*0.0,0.0};
    //Is really { {V[0]*0,0.0}, {0.0,0.0}, {0.0,0.0}, {1.0,0.0} };
    matrices[0](0,3)={1.0,0.0};

    matrices[1](0,1) = -t;//sig=0, sig1=1
    matrices[2](0,2) = -conj(t);//sig=1, sig1=0
    matrices[3](0,0) = {V[0],0.0};//sig=1, sig1=1
    matrices[3](0,3) = {1.0,0.0};//sig=1, sig1=1

    for (uint j = 1; j+1 < L; ++j)
    {

        matrices[0+dTwo*j](0,0)={1.0,0.0};
        matrices[0+dTwo*j](3,0)={V[j]*0,0.0};matrices[0+dTwo*j](3,3)={1.0,0.0};

        matrices[1+dTwo*j](2,0)={1.0,0.0};
                                  matrices[1+dTwo*j](3,1)=-t;

        matrices[2+dTwo*j](1,0)={1.0,0.0};
                                  matrices[2+dTwo*j](3,2)=-conj(t);


        matrices[3+dTwo*j](0,0)={1.0,0.0};
        matrices[3+dTwo*j](3,0)={V[j],0.0};matrices[3+dTwo*j](3,3)={1.0,0.0};
    }


    column_lookup[0]=0;
    column_lookup[1]=0;
    column_lookup[2]=0;


    matrices[L*dTwo-4](0,0)={1.0,0.0};
    matrices[L*dTwo-4](3,0)={V[L-1]*0,0.0};
    matrices[L*dTwo-3](2,0)={1.0,0.0};
    matrices[L*dTwo-2](1,0)={1.0,0.0};
    matrices[L*dTwo-1](0,0)={1.0,0.0};
    matrices[L*dTwo-1](3,0)={V[L-1],0.0};
}

mpo::mpo(
        uint64_t Lx,
        uint64_t Ly,
        //uint64_t N_particles,
        //uint64_t N_anyons,
        vector<double>& V,
        vector<complex<double> >& tx,
        vector<complex<double> >& ty
        )
{
    L=Lx*Ly;
    d=2;
    dTwo=d*d;

    Dw=4+2*(Ly);//Number of Jumps

    column_lookup=vector<uint>(Dw-1,0);//Not really worth anything here

    matrices=vector<sp_cx_mat>(L*dTwo,sp_cx_mat(Dw,Dw));//There are L*(d**2) matrices

    matrices[0]=sp_cx_mat(1,Dw);
    matrices[1]=sp_cx_mat(1,Dw);
    matrices[2]=sp_cx_mat(1,Dw);
    matrices[3]=sp_cx_mat(1,Dw);

    matrices[L*dTwo-4]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-3]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-2]=sp_cx_mat(Dw,1);
    matrices[L*dTwo-1]=sp_cx_mat(Dw,1);






    //sig=0, sig1=0
    matrices[0](0,0)={V[0]*0.0,0.0};
    //Is really { {V[0]*0,0.0}, {0.0,0.0}, {0.0,0.0}, {1.0,0.0} };
    matrices[0](0,Dw-1)={1.0,0.0};

    matrices[1](0,Dw-3)  = -ty[0];//sig=0, sig1=1
    matrices[1](0,1+Ly-1)= -tx[0];//sig=0, sig1=1

    matrices[2](0,Dw-2) = -conj(ty[0]);//sig=1, sig1=0
    matrices[2](0,1+2*Ly-1) = -conj(tx[0]);//sig=1, sig1=0

    matrices[3](0,0) = {V[0],0.0};//sig=1, sig1=1
    matrices[3](0,Dw-1) = {1.0,0.0};//sig=1, sig1=1


    //The entries from 0,0
    column_lookup[0]=0;
    column_lookup[Ly+1]=0;
    column_lookup[Dw-2]=0;

    for (uint i = 0; i<Ly-1; ++i)
    {
        column_lookup[i+2]=i+1;
        column_lookup[i+2+Ly]=i+1+Ly;
    }


    column_lookup[1]=0;
    column_lookup[1+Ly]=0;
    column_lookup[Dw-2]=0;
    column_lookup[Dw-3]=0;

    for (uint j = 1; j+1 < L; ++j)
    {
        //sig=0, sig1=0
        matrices[0+dTwo*j](0,0)   ={1.0   ,0.0};
        matrices[0+dTwo*j](Dw-1,0)={V[j]*0,0.0};matrices[0+dTwo*j](Dw-1,Dw-1)={1.0,0.0};


        for (uint i = 0; i<Ly-1; ++i)
        {
            matrices[0+dTwo*j](i+2   ,(i+1))    = {1.0,0.0};//sig=0, sig1=0
            matrices[0+dTwo*j](i+2+Ly,(i+1+Ly)) = {1.0,0.0};//sig=0, sig1=0
        }


        matrices[1+dTwo*j](Dw-2,0)={1.0,0.0};
        matrices[1+dTwo*j](1+Ly,0)={1.0,0.0};
                            matrices[1+dTwo*j](Dw-1,Dw-3)=-ty[j]; matrices[1+dTwo*j](Dw-1,1+Ly-1)=-tx[j];

        matrices[2+dTwo*j](1,0)   ={1.0,0.0};
        matrices[2+dTwo*j](Dw-3,0)={1.0,0.0};
                            matrices[2+dTwo*j](Dw-1,Dw-2)=-conj(ty[j]); matrices[2+dTwo*j](Dw-1,1+2*Ly-1)=-conj(tx[j]);

        matrices[3+dTwo*j](0,0)={1.0,0.0};
        matrices[3+dTwo*j](Dw-1,0)   ={V[j],0.0};matrices[3+dTwo*j](Dw-1,Dw-1)={1.0,0.0};

        for (uint i = 0; i<Ly-1; ++i)
        {
            matrices[3+dTwo*j](i+2,(i+1))      = {1.0,0.0};//sig=0, sig1=0
            matrices[3+dTwo*j](i+2+Ly,(i+1+Ly))= {1.0,0.0};//sig=0, sig1=0
        }
    }





    matrices[(L-1)*dTwo](0,0)        ={1.0,0.0};
    matrices[(L-1)*dTwo](Dw-1,0)     ={V[L-1]*0,0.0};

    matrices[(L-1)*dTwo+1](Dw-2,0)   ={1.0,0.0};
    matrices[(L-1)*dTwo+1](1+Ly,0)   ={1.0,0.0};

    matrices[(L-1)*dTwo+2](1,0)      ={1.0,0.0};
    matrices[(L-1)*dTwo+2](Dw-3,0)   ={1.0,0.0};

    matrices[(L-1)*dTwo+3](0,0)      ={1.0,0.0};
    matrices[(L-1)*dTwo+3](Dw-1,0)   ={V[L-1],0.0};
}




complex<double> mpo::get_expectation_value(const mps& mps1, uint j, bool verbose)
{

    if (mps1.can_center==-2)
        throw std::runtime_error("No canonization, can not proceed");


    //Now verify the expectation value
    uint n_c = mps1.can_center;

    //Pretend normalization is somewhere within the range were we can work with it
    if (mps1.can_center==L)
        n_c = L-1;
    if (mps1.can_center==-1)
        n_c = 0;
    //cout<<"Set norm center "<<n_c<<':'<<mps1.can_center<<endl;

    //Start from the left, move towards the center and get all the L's, then start from the right and move towards center to get all R's
    vector<cx_cube> tensors(L+1);
    tensors[0]=cx_cube(1,1,1,fill::ones);
    tensors[L]=cx_cube(1,1,1,fill::ones);


    //If this is the n operator
    if (is_n && j < L)
    {
        //Basically just copy in 1 (the identity) at a particular point in the data, and copy it away afterwards
        if (j==0)
            matrices[3](0,0) = {1.0,0.0};
        else if (j+1==L)
            matrices[L*dTwo-1](1,0)={1.0,0.0};
        else
            matrices[3+dTwo*j](1,0)={1.0,0.0};
    }


    auto start = chrono::steady_clock::now();

    if (n_c>=1)
    {
        for (uint l = 1; l<=n_c; ++l)
        {
            tensors[l]=getL_recursive(mps1,l, tensors[l-1]);

        }
    }


    if (n_c<L)
    {
        for (uint l = L-1; l>n_c && l<L; --l)
        {
            tensors[l]=getR_recursive(mps1,l, tensors[l+1]);
        }
    }

    auto endLR = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(endLR - start);

    // To get the value of duration use the count()
    // member function on the duration object
    if (verbose)
        cout <<"Time to LR "<< duration.count() <<" μs"<<endl;


/*


    cout<<"L GOTTEN "<<endl;

    for (uint l = 0; l<=n_c; ++l)
    {
        cout<<l<<endl;
        //Previous and current m,n of the mps (a) and MPO (b)
        uint mAl  = mps1.D_list[l];
        uint mWl  = l==0 ? 1 : Dw;

        for (uint a0 = 0; a0 < mAl; ++a0)
        {
            for (uint a1 = 0; a1 < mAl; ++a1)
            {

                cout<<a0<<" "<<a1<<endl;

                for (uint b = 0; b < mWl; ++b)
                {

                    cout<<tensors[l](a0,a1,b)<<' ';
                }
                cout<<endl;
            }
        }
    }

    cout<<"R GOTTEN "<<endl;

    for (uint l = n_c+1; l<=L; ++l)
    {
        cout<<l<<endl;
        //Previous and current m,n of the mps (a) and MPO (b)
        uint mAl  = mps1.D_list[l];
        uint mWl  = l==L ? 1 : Dw;

        for (uint a0 = 0; a0 < mAl; ++a0)
        {
            for (uint a1 = 0; a1 < mAl; ++a1)
            {

                cout<<a0<<" "<<a1<<endl;
                for (uint b = 0; b < mWl; ++b)
                {
                    cout<<" "<<b<<' '<<flush;


                    cout<<tensors[l](a0,a1,b)<<' ';
                }
                cout<<endl;
            }
        }
    }
*/


    uint   a_min=0;
    uint   a_max=0;
    uint ap1_min=0;
    uint ap1_max=0;

    vector<submatrix_index> submatrix_indices;

    for (uint sig = 0 ; sig<d; ++sig)
        for (uint k = 0; mps1.sub_matrix_range(n_c, sig, k, a_min, a_max, ap1_min, ap1_max) ; ++k)
            for (uint a = a_min; a<a_max;++a)
                for (uint ap1 = ap1_min; ap1<ap1_max;++ap1)
        {
            submatrix_indices.push_back(submatrix_index(sig, k, a, ap1,a_min,ap1_min));
        }


    //Sum over all tensors, the current MPS and current MPO matrices to get expectation value

    //Size of current MPS
    uint mW = n_c==0? 1 : Dw;
    uint nW = n_c+1==L? 1 : Dw;

    complex<double> expectation_value ={0.0,0.0};

    //We could also get this matrix, and this vector
    //sp_cx_mat grand_H(submatrix_indices.size(),submatrix_indices.size());
    //cx_vec M_vector(submatrix_indices.size());

    //Skip anything if the MPO matrix is 0 (mostly empty)
    for (uint I = 0; I<submatrix_indices.size(); ++I)
    {
       // M_vector(I)=mps1.get(submatrix_indices[I].sig,n_c, submatrix_indices[I].k,submatrix_indices[I].a-submatrix_indices[I].a_min,submatrix_indices[I].ap1-submatrix_indices[I].ap1_min);
        for (uint J = 0; J<submatrix_indices.size(); ++J)
        {
            for (uint bp1 = 0; bp1<nW; ++bp1)//One of the coordinates to write to
                for (uint b = 0; b<mW; ++b)
                    if (!(approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).real(),0.0) && approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).imag(),0.0)) )
                    {
                        //Have to cast this to a complex for some stupid reason, as it otherwise gets interpreted as a matrix
                        complex<double> NOT_MATRICE =matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1);
                        //grand_H(I,J)+=NOT_MATRICE*tensors[n_c](submatrix_indices[I].a,submatrix_indices[J].a,b)*tensors[n_c+1](submatrix_indices[I].ap1,submatrix_indices[J].ap1,bp1);
                           expectation_value  +=NOT_MATRICE *mps1.get(submatrix_indices[J].sig,n_c, submatrix_indices[J].k,submatrix_indices[J].a-submatrix_indices[J].a_min,submatrix_indices[J].ap1-submatrix_indices[J].ap1_min)*tensors[n_c](submatrix_indices[I].a,submatrix_indices[J].a,b)*tensors[n_c+1](submatrix_indices[I].ap1,submatrix_indices[J].ap1,bp1)*conj(mps1.get(submatrix_indices[I].sig,n_c, submatrix_indices[I].k,submatrix_indices[I].a-submatrix_indices[I].a_min,submatrix_indices[I].ap1-submatrix_indices[I].ap1_min));
                    }

        }
    }


    //We could have gotten the result like this as well, but that is slower
    /*
    cx_vec prod = grand_H*M_vector;

    complex<double> res = dot(M_vector.t(),prod);
    */

    auto endSum = chrono::steady_clock::now();

    duration = chrono::duration_cast<chrono::microseconds>(endSum -endLR);


    //If this is the n operator, undo the changes we just did
    if (is_n && j < L)
    {
        //Basically just copy in 1 (the identity) at a particular point in the data, and copy it away afterwards
        if (j==0)
            matrices[3](0,0) = {0.0,0.0};
        else if (j+1==L)
            matrices[L*dTwo-1](1,0)={0.0,0.0};
        else
            matrices[3+dTwo*j](1,0)={0.0,0.0};
    }



    // To get the value of duration use the count()
    // member function on the duration object
    if (verbose)
        cout <<"Time to sum"<< duration.count() <<" μs" << endl;
    return expectation_value;
}

//Get the L[l] tensor, calculated given the tensor L[l-1]. L[l] has indices a0_l, a1_l and b_l, where  a_l are the indices of the left matrices of the MPO
cx_cube mpo::getL_recursive(const mps& mps1, uint j_plus1,const cx_cube& L_prev) const
{
    if (j_plus1==0)
        return cx_cube(1,1,1,fill::ones);
    //j is specificall

    uint j = j_plus1-1;

    uint mAj  = mps1.D_list[j];
    uint mWj  = j==0 ? 1 : Dw;


    uint nAj  = mps1.D_list[j_plus1];
    uint nWj  = j_plus1==L ? 1 : Dw;


    cx_cube L_tensor(nAj,nAj,nWj,fill::zeros);

    //We want to quickly eliminate as many Loop as possible, and all matrices in the multiplications are mostly 0 anyway

    //From my testing, this appears to be the fastest way to nest these loops, and the best way to reuse as many calculations as we can get away with the basic idea is to discount as much as we can as early as possible, remember most of everything is very very sparse

    //Skip anything if the MPO matrix is 0 (mostly empty)
    for (uint sig0=0; sig0<d; ++sig0)
    {
        uint a0j_min=0;
        uint a0j_max=mAj;
        uint a0jp1_min=0;
        uint a0jp1_max=nAj;

        uint a1j_min=0;
        uint a1j_max=mAj;
        uint a1jp1_min=0;
        uint a1jp1_max=nAj;

        //This function is somewhat expensive, so we will repeat it as few times as need be, we need sig0 (but it is absolutely worth it, because it picks out the submatrices alone; looping only inside them gives a x10 improvement; I have considered pre-calculating these things, but that does not seem worth it)
        for (uint k0 = 0; mps1.sub_matrix_range(j, sig0, k0, a0j_min, a0j_max, a0jp1_min, a0jp1_max) ; ++k0)
        {
            for (uint sig1=0; sig1<d; ++sig1)
                for (uint k1 = 0; mps1.sub_matrix_range(j, sig1, k1, a1j_min, a1j_max, a1jp1_min, a1jp1_max) ; ++k1)
                {
                    for (uint bjp1 = 0; bjp1<nWj; ++bjp1)//One of the coordinates to write to
                        for (uint bj = 0; bj<mWj; ++bj)
                            if (!(approx(matrices[sig0*d+sig1+dTwo*j](bj,bjp1).real(),0.0) && approx(matrices[sig0*d+sig1+dTwo*j](bj,bjp1).imag(),0.0)) )
                            {
                                for (uint a0jp1 = a0jp1_min; a0jp1<a0jp1_max; ++a0jp1)
                                    for (uint a1jp1 = a1jp1_min; a1jp1<a1jp1_max; ++a1jp1)
                                    {
                                        complex<double> MLM_sum={0.0,0.0};
                                        //TEMP, replace with loop through submatrices
                                        for (uint a0j= a0j_min; a0j<a0j_max; ++a0j)//TEMP, check if 'tis 0
                                        {
                                            complex<double> sum_over_allML={0.0,0.0};
                                            for (uint a1j= a1j_min; a1j<a1j_max; ++a1j)
                                            {
                                                sum_over_allML+=mps1.get(sig1,j, k1,a1j-a1j_min,a1jp1-a1jp1_min)*L_prev(a0j,a1j,bj);
                                            }
                                            MLM_sum += conj(mps1.get(sig0,j, k0,a0j-a0j_min,a0jp1-a0jp1_min))*sum_over_allML;
                                        }
                                        L_tensor(a0jp1,a1jp1,bjp1)+=MLM_sum*matrices[sig0*d+sig1+dTwo*j](bj,bjp1);
                                    }
                            }
                }
        }
    }
    return L_tensor;
}
//Same for R[l]
cx_cube mpo::getR_recursive(const mps& mps1, uint j,const cx_cube& R_next) const
{
    if (j==L)
        return cx_cube(1,1,1,fill::ones);
    //j is specificall


    uint j_plus1 = j+1;


    uint mAj  = mps1.D_list[j];
    uint mWj  = j==0 ? 1 : Dw;


    uint nAj  = mps1.D_list[j_plus1];
    uint nWj  = j_plus1==L ? 1 : Dw;

    cx_cube R_tensor(mAj,mAj,mWj,fill::zeros);


    //From my testing, this appears to be the fastest way to nest these loops, and the best way to reuse as many calculations as we can get away with the basic idea is to discount as much as we can as early as possible, remember most of everything is very very sparse
    for (uint sig0=0; sig0<d; ++sig0)
    {
        uint a0j_min=0;
        uint a0j_max=mAj;
        uint a0jp1_min=0;
        uint a0jp1_max=nAj;

        uint a1j_min=0;
        uint a1j_max=mAj;
        uint a1jp1_min=0;
        uint a1jp1_max=nAj;

        //This function is somewhat expensive, so we will repeat it as few times as need be, we need sig0 (but it is absolutely worth it, because it picks out the submatrices alone; looping only inside them gives a x10 improvement; I have considered pre-calculating these things, but that does not seem worth it)
        for (uint k0 = 0; mps1.sub_matrix_range(j, sig0, k0, a0j_min, a0j_max, a0jp1_min, a0jp1_max) ; ++k0)
        {
            for (uint sig1=0; sig1<d; ++sig1)
                for (uint k1 = 0; mps1.sub_matrix_range(j, sig1, k1, a1j_min, a1j_max, a1jp1_min, a1jp1_max) ; ++k1)
                {
                    for (uint bjp1 = 0; bjp1<nWj; ++bjp1)//One of the coordinates to write to
                        for (uint bj = 0; bj<mWj; ++bj)
                            if (!(approx(matrices[sig0*d+sig1+dTwo*j](bj,bjp1).real(),0.0) && approx(matrices[sig0*d+sig1+dTwo*j](bj,bjp1).imag(),0.0)) )
                            {
                                for (uint a0j= a0j_min; a0j<a0j_max; ++a0j)//TEMP, check if 'tis 0
                                    for (uint a1j= a1j_min; a1j<a1j_max; ++a1j)
                                    {
                                        complex<double> MLM_sum={0.0,0.0};
                                        //TEMP, replace with loop through submatrices
                                        for (uint a0jp1 = a0jp1_min; a0jp1<a0jp1_max; ++a0jp1)
                                        {
                                            complex<double> sum_over_allML={0.0,0.0};
                                            for (uint a1jp1 = a1jp1_min; a1jp1<a1jp1_max; ++a1jp1)
                                            {
                                                sum_over_allML+=mps1.get(sig1,j, k1,a1j-a1j_min,a1jp1-a1jp1_min)*R_next(a0jp1,a1jp1,bjp1);
                                            }
                                            MLM_sum += conj(mps1.get(sig0,j, k0,a0j-a0j_min,a0jp1-a0jp1_min))*sum_over_allML;
                                        }
                                        R_tensor(a0j,a1j,bj)+=MLM_sum*matrices[sig0*d+sig1+dTwo*j](bj,bjp1);
                                    }
                            }
                }
        }
    }

    return R_tensor;

}


void mpo::DEBUG_PRINT_ALL() const
{
    for (uint j = 0; j<L; ++j)
    {
        cout<<" j = "<<j<<endl;
        for (uint sig0 = 0; sig0<d; ++sig0)
            for (uint sig1 = 0; sig1<d; ++sig1)
            {
                cout<<sig0<<' '<<sig1<<endl;
                for (uint i = 0; i < (j==0 ? 1 : Dw); ++i)
                {
                    for (uint k = 0; k < (j+1==L ? 1 : Dw); ++k)
                        cout<<matrices[sig0*d+sig1+dTwo*j](i,k)<<' ';
                    cout<<endl;
                }
            }

    }
}


mps mpo::find_ground(const mps& that, double& energy,uint max_sweeps)
{

    //cout<<"Start Find ground"<<endl;
    //Verify length and d matches
    if (that.get_d()!=d || that.get_L()!=L)
        throw std::runtime_error("matrices do not have same length or dimension");

    //Copy, but extend the mps dimenusion

    mps new_mps = that;//Default copy constructor is fine (I think)

    if (new_mps.can_center !=0)
    {
        //cout<<"Right-canonizing..."<<endl;
        while(new_mps.step_left_QR()){}
        //cout<<"Done"<<endl;
    }


    uint n_c = new_mps.can_center;
    uint sweep = 0;//Which sweep are we on

    //For looking at convergence, energy at each step
    vector<double> Energies;
    Energies.reserve(max_sweeps*2*(L-1));

    double E_threshold=1e-4;//If less than this change for a single half-sweep, we break early
    double E_sweepstart=0;

    uint step_id=0;

    bool go_right=true;


    //cout<<"Recursively creating All tensors "<<endl;


    //Start from the left, move towards the center and get all the L's, then start from the right and move towards center to get all R's
    vector<cx_cube> tensors(L+1);
    tensors[0]=cx_cube(1,1,1,fill::ones);
    tensors[L]=cx_cube(1,1,1,fill::ones);

    auto start = chrono::steady_clock::now();

    if (n_c>=1)
    {
        for (uint l = 1; l<=n_c; ++l)
        {
            tensors[l]=getL_recursive(new_mps,l, tensors[l-1]);

        }
    }


    if (n_c<L)
    {
        for (uint l = L-1; l>n_c && l<L; --l)
        {
            tensors[l]=getR_recursive(new_mps,l, tensors[l+1]);
        }
    }

    auto endLR = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(endLR - start);


    vector<uint> process_counts;




    //may return 0 when not able to detect
    const auto processor_count = std::thread::hardware_concurrency();

    cout<<"Using "<<processor_count<<" threads "<< endl;



    while(n_c<L && sweep<max_sweeps)
    {


        auto start = chrono::steady_clock::now();




        auto start_sweep = chrono::steady_clock::now();
        //cout<<"Sweep "<<sweep<<" Site "<<n_c<<"/"<<(L-1)<<" going "<<step_id<<endl;;
/*
        {
            cout<<"Sweep "<<sweep<<" Site "<<n_c<<"/"<<(L-1)<<" going";
            if (go_right)
                cout<<" right"<<endl;
            else
                cout<<" left"<<endl;
        }*/
        //Previous and current m,n of the mps (M) and MPO (W)

        uint   a_min=0;
        uint   a_max=0;
        uint ap1_min=0;
        uint ap1_max=0;

        vector<submatrix_index> submatrix_indices;

        for (uint sig = 0 ; sig<d; ++sig)
            for (uint k = 0; new_mps.sub_matrix_range(n_c, sig, k, a_min, a_max, ap1_min, ap1_max) ; ++k)
                for (uint a = a_min; a<a_max;++a)
                    for (uint ap1 = ap1_min; ap1<ap1_max;++ap1)
                        submatrix_indices.push_back(submatrix_index(sig, k, a, ap1,a_min,ap1_min));
       uint reduced_size=submatrix_indices.size();



        sp_cx_mat grand_H(reduced_size,reduced_size);

        //Size of current MPS
        uint mW = n_c==0? 1 : Dw;
        uint nW = n_c+1==L? 1 : Dw;

   //cout<<" REDUCED SIZE GOTTEN "<<reduced_size<<" MPO SIZE IS GOT "<<nW<<" "<<mW<<endl;

        auto start_H = chrono::steady_clock::now();

/* 1 thread version, better for smaller systems
        //Skip anything if the MPO matrix is 0 (mostly empty)

        for (uint I = 0; I<reduced_size; ++I)
        {
            for (uint J = 0; J<reduced_size; ++J)
            {
                uint b=0;
                for (b = 0; b<mW-1; ++b)
                {
                    uint bp1 = column_lookup[b];

                    if (bp1<nW)
                    {

                        if (!(approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).real(),0.0) && approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).imag(),0.0)) )
                        {
                            //Have to cast this to a complex for some stupid reason, as it otherwise gets interpreted as a matrix
                            complex<double> NOT_MATRICE =matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1);
                            grand_H(I,J)+=NOT_MATRICE*tensors[n_c](submatrix_indices[I].a,submatrix_indices[J].a,b)*tensors[n_c+1](submatrix_indices[I].ap1,submatrix_indices[J].ap1,bp1);

                        }
                    }
                }
                b =mW-1;
                for (uint bp1 = 0; bp1<nW; ++bp1)
                {

                    if (!(approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).real(),0.0) && approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).imag(),0.0)) )
                    {
                        //Have to cast this to a complex for some stupid reason, as it otherwise gets interpreted as a matrix
                        complex<double> NOT_MATRICE =matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1);
                        grand_H(I,J)+=NOT_MATRICE*tensors[n_c](submatrix_indices[I].a,submatrix_indices[J].a,b)*tensors[n_c+1](submatrix_indices[I].ap1,submatrix_indices[J].ap1,bp1);

                    }
                }
            }
        }
*/
        uint threads = processor_count;
        uint IsPer_thread = reduced_size/threads;

        vector <sp_cx_mat> not_so_grand_H0(threads,sp_cx_mat(reduced_size,reduced_size));

        vector<std::thread> Threads;

        for (uint i = 0; i <threads; ++i)
        {
             uint I_min = i*IsPer_thread;
             uint I_max = (i == threads-1 ? reduced_size : (i+1)*IsPer_thread );
            //make_H_thread(I_min , I_max, reduced_size,mW, nW,n_c, submatrix_indices, tensors[n_c], tensors[n_c+1], &grand_H, d, dTwo, column_lookup, matrices);

            Threads.push_back(std::thread(make_H_thread, I_min , I_max, reduced_size,mW, nW,n_c, submatrix_indices, tensors[n_c], tensors[n_c+1], &not_so_grand_H0[i], d, dTwo, column_lookup, matrices));
        }

        for (uint i = 0; i <threads; ++i)
        {
            Threads[i].join();
            grand_H+=not_so_grand_H0[i];
        }



        auto has_H = chrono::steady_clock::now();
        //Solve eigenvalue

        //cout<<" sum cost "<<chrono::duration_cast<chrono::microseconds>(has_H - start_H ).count()<<" other "<<chrono::duration_cast<chrono::microseconds>(start_H - start_sweep ).count()<<endl;


        //eigs_gen = [eig]envalues of [s]parse [gen]eral (i.e. need not be real) matrix

        cx_vec eigval(1);
        cx_mat eigvec(reduced_size,1);

        uint smallest_eigenvalue = 0;//If we use the sparse eigensolver, the smallest eigenvalue is the only found, only in a single stupid special case might that not be so

        if (reduced_size==2)//This is the smallest we can get... and it is terrible! because the solver breaks
        {
            //Oh well, bring out the dense solver, a 2 by 2 isn't that big anyway
            cx_mat tmp(grand_H);

            eig_gen( eigval, eigvec, tmp);

            if (eigval[0].real()>eigval[1].real())
                smallest_eigenvalue=1;

        }
        else
        {
            eigs_gen( eigval, eigvec, grand_H, 1, "sr");//The lowest eigenstate
        }


        auto has_solved = chrono::steady_clock::now();

        energy=eigval(smallest_eigenvalue).real();
        //cout<<" --- Energy "<<n_c<<' '<<energy<<" Using size "<<reduced_size<<endl;


        Energies.push_back(energy);


        if (step_id==0)//This is the first energy
            E_sweepstart=energy;
        ++step_id;

        //Now read off the new matrix

        for (uint I = 0; I < reduced_size; ++I)
        {
            new_mps.set(submatrix_indices[I].sig,n_c, submatrix_indices[I].k,submatrix_indices[I].a-submatrix_indices[I].a_min,submatrix_indices[I].ap1-submatrix_indices[I].ap1_min,eigvec(I,smallest_eigenvalue));

        }


        //DEBUG
        //cout<<"CHECK EXPECTATION VALUE"<<endl;
        //cout<<get_expectation_value(new_mps)<<endl;

        auto read_newmatrix= chrono::steady_clock::now();
        //MORE DEBUG

        //To prepare for the next step, step either left or right
        if (go_right && n_c+1==L)
        {
            go_right =false;

            //cout<<"Sweep "<<sweep<<"right done, with energy "<<energy<<endl;

            if (step_id!=0)
                if (std::abs(energy-E_sweepstart) <E_threshold)
                {
                    cout<<"Energy converged throughout sweep, breaking in middle"<<endl;
                    break;
                }


            E_sweepstart=energy;
        }
        else if (!go_right  && n_c==0)
        {
            //cout<<"Sweep left "<<sweep<<" done, with energy "<<energy<<endl;
            if (step_id!=0)
                if (std::abs(energy-E_sweepstart) <E_threshold)
                {
                    cout<<"Energy converged throughout sweep, breaking at end"<<endl;
                    break;
                }
            go_right =true;
            ++sweep;
            E_sweepstart=energy;
        }

        process_counts.push_back(chrono::duration_cast<chrono::microseconds>(has_H - start_sweep).count());
        process_counts.push_back(chrono::duration_cast<chrono::microseconds>(has_solved  - has_H ).count());
        if (go_right)
        {
           // cout<<"Steping Right "<<energy<<endl;
            new_mps.step_right_QR();//Step right moves norm center right, so this is a make left
            auto QR = chrono::steady_clock::now();
            //cout<<"Update Tensor"<<endl;

//DO NOT DO THIS
            if (n_c+1<L)
                tensors[n_c+1]=getL_recursive(new_mps,n_c+1, tensors[n_c]);


            ++n_c;

            auto getL = chrono::steady_clock::now();

            process_counts.push_back(chrono::duration_cast<chrono::microseconds>(QR - has_solved ).count());
            process_counts.push_back(chrono::duration_cast<chrono::microseconds>(getL - QR ).count());


        }
        else
        {
        //    cout<<"Steping Left "<<energy<<endl;
            new_mps.step_left_QR();//And vica versa

            auto QR = chrono::steady_clock::now();
            //cout<<"Update Tensor"<<endl;

            tensors[n_c]=getR_recursive(new_mps,n_c, tensors[n_c+1]);//The final is alway [[[1]]]
        --n_c;


            auto getR = chrono::steady_clock::now();



            process_counts.push_back(chrono::duration_cast<chrono::microseconds>(QR - has_solved ).count());
            process_counts.push_back(chrono::duration_cast<chrono::microseconds>(getR - QR ).count());
        }

    }

    cout<<"Performed "<<sweep<<" sweeps (max "<<max_sweeps<<")"<<endl;

    cout <<"Time to LR "<< duration.count()  <<" μs" << endl;
    uint   H_time_total=0;
    uint eig_time_total=0;
    uint  QR_time_total=0;
    uint ten_time_total=0;

    uint time_total=0;
    uint time2_total=0;
    uint steps=0;
    for (uint i = 0; i < process_counts.size(); i+=4)
    {
        ++steps;
        //cout<<"Step "<<i/4<<endl;
        //cout<<"Setting up H matrix:"<<process_counts[i+0]<<"μs"<<endl ;
        H_time_total+=process_counts[i+0];
        //cout<<"Get eigenvalues    :"<<process_counts[i+1]<<"μs"<<endl ;
        eig_time_total+=process_counts[i+1];
        //cout<<"Use QR canonization:"<<process_counts[i+2]<<"μs"<<endl ;
        QR_time_total+=process_counts[i+2];
        //cout<<"Update tensors     :"<<process_counts[i+3]<<"μs"<<endl ;
        ten_time_total+=process_counts[i+3];

        uint t = process_counts[i+0]+process_counts[i+1]+process_counts[i+2]+process_counts[i+3];
        time_total+=t;
    }

    uint avg = time_total/steps;

    for (uint i = 0; i < process_counts.size(); i+=4)
    {

        uint t = process_counts[i+0]+process_counts[i+1]+process_counts[i+2]+process_counts[i+3];
        time2_total+=int(t-avg)*int(t-avg);
    }






    cout<<"AVERAGED time cost over "<<steps<<" steps "<<endl;
    cout<<"Setting up H matrix:"<<H_time_total/steps<<"μs"<<endl ;
    cout<<"Get eigenvalues    :"<<eig_time_total/steps<<"μs"<<endl ;
    cout<<"Use QR canonization:"<<QR_time_total/steps<<"μs"<<endl ;
    cout<<"Update tensors     :"<<ten_time_total/steps<<"μs"<<endl ;

    cout<<time2_total<<endl;

    cout<<"Total "<<avg<<" +- "<<sqrt(time2_total/steps)<<endl;


    auto end = chrono::steady_clock::now();

/*
    cout<<"Energies: "<<endl;
    for (double& energies : Energies)
    {
        cout<<energies <<endl;
    }*/
    cout<<"Total time: "<<chrono::duration_cast<chrono::microseconds>(end - start).count()<<" μs"<<endl;

    return new_mps;
}

void mpo::set_potential(const vector<double>& V)
{

    matrices[3](0,0)={V.size()>0 ? V[0] : 0.0,0.0};

    for (uint j = 1; j<L ; ++j)
        matrices[3+dTwo*j](Dw-1,0)={V.size()>0 ? V[j] : 0.0,0.0};

}

void mpo::make_H_thread(uint I_min, uint I_max, uint reduced_size,uint mW, uint nW,uint n_c,const vector<submatrix_index>& submatrix_indices, const cx_cube& tensorL, const cx_cube& tensorR, sp_cx_mat* grand_H,uint d, uint dTwo, const vector<uint>& column_lookup, const vector<sp_cx_mat>& matrices)
{
    for (uint I = I_min; I<I_max; ++I)
    {
        for (uint J = 0; J<reduced_size; ++J)
        {
            uint b=0;
            for (b = 0; b<mW-1; ++b)
            {
                uint bp1 = column_lookup[b];

                if (bp1<nW)
                {

                    if (!(approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).real(),0.0) && approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).imag(),0.0)) )
                    {
                        //Have to cast this to a complex for some stupid reason, as it otherwise gets interpreted as a matrix
                        complex<double> NOT_MATRICE =matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1);
                        (*grand_H)(I,J)+=NOT_MATRICE*tensorL(submatrix_indices[I].a,submatrix_indices[J].a,b)*tensorR(submatrix_indices[I].ap1,submatrix_indices[J].ap1,bp1);

                    }
                }
            }
            b =mW-1;
            for (uint bp1 = 0; bp1<nW; ++bp1)
            {

                if (!(approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).real(),0.0) && approx(matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1).imag(),0.0)) )
                {
                    //Have to cast this to a complex for some stupid reason, as it otherwise gets interpreted as a matrix
                    complex<double> NOT_MATRICE =matrices[submatrix_indices[I].sig*d+submatrix_indices[J].sig+dTwo*n_c](b,bp1);
                    (*grand_H)(I,J)+=NOT_MATRICE*tensorL(submatrix_indices[I].a,submatrix_indices[J].a,b)*tensorR(submatrix_indices[I].ap1,submatrix_indices[J].ap1,bp1);

                }
            }
        }
    }
}
