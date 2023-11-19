#include "mps.hpp"
#include "approx.hpp"



uint64_t bin_coef(uint64_t n, uint64_t m);

//Uniform starting guess
mps::mps(uint64_t _L,uint _N_particles, uint local_D, bool left)
{

    L =_L;
    d = 2;

    N_particles=_N_particles;
    can_center=-1;

    //These all include a final "faux" matrix with height 1, to make the width of the final thing 1
    D_list             = vector<uint >(L+1,1);
    subdimension_height = vector<vector<uint> >(L+1,vector<uint>(1,1));

    row_start=vector<vector<uint> >(L+1,vector<uint>(1,0));

    n_subrows=vector<uint >(L+1,1);
    sub_matrices       = vector< vector<cx_mat> >(d*L);
    top_submatrix_location=vector<pair<uint,uint> >(d*L);
    firstsubcol_totalparticles=vector< uint > (L);

    uint max_subrow = 0;
    uint min_subrow = 0;

    for (uint j = 0; j < L; ++j)
    {

        //These are the bounds in the true faux matrix where our blocks exist, for the first and the last it is strictly speaking not true, but the test will still work

        uint max_subcolumn = (j+1>N_particles ? N_particles : j+1 );
        uint min_subcolumn = j+N_particles<L ?  0 : j*(d-1)-(d-1)*(L-1)+N_particles;

        firstsubcol_totalparticles[j]=min_subcolumn ;

        n_subrows[j+1]= max_subcolumn-min_subcolumn+1;

        subdimension_height[j+1]=vector<uint>(n_subrows[j+1]);
        row_start[j+1] =vector<uint>(n_subrows[j+1]);

        for (uint k = 0; k<n_subrows[j+1]; ++k)
        {

            if (j+1==L)
            {
                subdimension_height[j+1][k]=1;
                row_start[j+1][k]          =k;
            }
            else
            {
                subdimension_height[j+1][k]=local_D;
                row_start[j+1][k]          =local_D*k;
            }
        }
        if (j+1==L)
            D_list[j+1]=1;
        else
            D_list[j+1]=n_subrows[j+1]*local_D;

        for (uint sig = 0; sig<d; ++sig)
        {
            //If there was no cut-off, this would be out top submatrix:
            //top_submatrix_location[d*j+sig]=make_pair(0,sig);
            //We pass the minimum subrow first
            if (min_subrow+sig<min_subcolumn)
            {
                top_submatrix_location[d*j+sig]=make_pair(min_subcolumn-min_subrow-sig,0);

            }
            else
            {
                top_submatrix_location[d*j+sig]=make_pair(0,-min_subcolumn+min_subrow+sig);

            }

            sub_matrices[j*d+sig]=vector<cx_mat>(std::min(1+max_subrow-min_subrow-top_submatrix_location[d*j+sig].first,1+max_subcolumn-min_subcolumn-top_submatrix_location[d*j+sig].second));

            double base = 1.0/sqrt(bin_coef(L,N_particles));
            for (uint k =0; k<sub_matrices[j*d+sig].size(); ++k)
            {
                sub_matrices[j*d+sig][k] = zeros<cx_mat>(subdimension_height[j][k],subdimension_height[j+1][k]);

                sub_matrices[j*d+sig][k](0,0) = (j==0 ? complex<double> {base,0.0} : complex<double>{1.0,0.0});
            }
        }

        max_subrow = max_subcolumn;
        min_subrow = min_subcolumn;
    }

    can_center=-2;
    if (!left)
    {
        //Start from the right and make everything right
        while (step_left_QR())
        {
        }

    }
    else
    {
        //Start from the left and make everything left
        while (step_right_QR())
        {
        }
    }

}

mps::mps(uint64_t _L,uint _N_particles,
const vector<pair<vector<bool> ,complex<double> >
    >& states,uint min_D, bool left)
{
    L =_L;
    d = 2;
    uint local_D = min_D>states.size()?  min_D : states.size() ;

    N_particles=_N_particles;
    can_center=-1;

    //These all include a final "faux" matrix with height 1, to make the width of the final thing 1
    D_list             = vector<uint >(L+1,1);
    subdimension_height = vector<vector<uint> >(L+1,vector<uint>(1,1));

    row_start=vector<vector<uint> >(L+1,vector<uint>(1,0));

    n_subrows=vector<uint >(L+1,1);
    sub_matrices       = vector< vector<cx_mat> >(d*L);
    top_submatrix_location=vector<pair<uint,uint> >(d*L);
    firstsubcol_totalparticles=vector< uint > (L);

    uint max_subrow = 0;
    uint min_subrow = 0;

    for (uint j = 0; j < L; ++j)
    {

        //These are the bounds in the true faux matrix where our blocks exist, for the first and the last it is strictly speaking not true, but the test will still work

        uint max_subcolumn = (j+1>N_particles ? N_particles : j+1 );
        uint min_subcolumn = j+N_particles<L ?  0 : j*(d-1)-(d-1)*(L-1)+N_particles;

        firstsubcol_totalparticles[j]=min_subcolumn ;

        n_subrows[j+1]= max_subcolumn-min_subcolumn+1;

        subdimension_height[j+1]=vector<uint>(n_subrows[j+1]);
        row_start[j+1] =vector<uint>(n_subrows[j+1]);

        for (uint k = 0; k<n_subrows[j+1]; ++k)
        {

            if (j+1==L)
            {
                subdimension_height[j+1][k]=1;
                row_start[j+1][k]          =k;
            }
            else
            {
                subdimension_height[j+1][k]=local_D;
                row_start[j+1][k]          =local_D*k;
            }
        }
        if (j+1==L)
            D_list[j+1]=1;
        else
            D_list[j+1]=n_subrows[j+1]*local_D;

        for (uint sig = 0; sig<d; ++sig)
        {
            //If there was no cut-off, this would be out top submatrix:
            //top_submatrix_location[d*j+sig]=make_pair(0,sig);
            //We pass the minimum subrow first
            if (min_subrow+sig<min_subcolumn)
            {
                top_submatrix_location[d*j+sig]=make_pair(min_subcolumn-min_subrow-sig,0);

            }
            else
            {
                top_submatrix_location[d*j+sig]=make_pair(0,-min_subcolumn+min_subrow+sig);

            }

            sub_matrices[j*d+sig]=vector<cx_mat>(std::min(1+max_subrow-min_subrow-top_submatrix_location[d*j+sig].first,1+max_subcolumn-min_subcolumn-top_submatrix_location[d*j+sig].second));

            for (uint k =0; k<sub_matrices[j*d+sig].size(); ++k)
            {
                sub_matrices[j*d+sig][k] = zeros<cx_mat>(subdimension_height[j][k],subdimension_height[j+1][k]);
                sub_matrices[j*d+sig][k] = zeros<cx_mat>(subdimension_height[j][k],subdimension_height[j+1][k]);

                for (uint l =0; l<states.size(); ++l)
                {
                    if ( (states[l].first[j] && sig==1 ) || (sig == 0 && !states[l].first[j]))
                    {
                        sub_matrices[j*d+sig][k]( j==0? 0:l, j==L-1? 0:l)= (j == 0 ? states[l].second : 1);
                    }
                }
            }
        }

        max_subrow = max_subcolumn;
        min_subrow = min_subcolumn;
    }

    can_center=-2;
    if (!left)
    {
        //Start from the right and make everything right
        while (step_left_QR())
        {
        }

    }
    else
    {
        //Start from the left and make everything left
        while (step_right_QR())
        {
        }
    }


}


mps::~mps()
{

}



void mps::DEBUG_PRINT_ALL() const
{
    std::cout << std::fixed << std::showpoint;
    //std::cout << std::setprecision(3);

    for (uint j = 0; j<L; ++j)
    {
        for (uint sig = 0; sig <d; ++sig)
            cout<<sub_matrices[j*d+sig].size()<<' ';
        cout<<endl;

        cout<<"Subrow number "<<n_subrows[j]<<endl;
        for (uint k = 0; k<n_subrows[j];++k)
        {

            for (uint m = 0; m<subdimension_height[j][k];++m)
            {

                for (uint sig = 0; sig <d; ++sig)
                {
                    for (uint l = 0; l<n_subrows[j+1];++l)
                    {
                        uint matrix_id = -1;

                        if (k>=top_submatrix_location[d*j+sig].first)
                        {
                            matrix_id =k-top_submatrix_location[d*j+sig].first;
                            if (top_submatrix_location[d*j+sig].second+matrix_id==l)
                            {


                            }
                            else
                                matrix_id =-1;
                        }

                        if (matrix_id == (uint)-1)
                            cout<<'|';
                        else
                            cout<<'[';
                        for (uint  n= 0; n<subdimension_height[j+1][l];++n)
                        {

                            if (matrix_id != (uint)-1)
                            {
                                cout<<sub_matrices[j*d+sig][matrix_id](m,n)<<' ';
                            }
                            else
                                cout<<(complex<double> {0.0,0.0})<<' ';

                        }

                        if (matrix_id == (uint)-1)
                            cout<<'|';
                        else
                            cout<<']';

                    }

                    cout<<"\t\t";
                }

                cout<<'\n';
            }
        }
        cout<<'\n';
        cout<<'\n';
    }

    cout <<flush;
}



complex<double> mps::get_coefficient(vector<bool> state) const
{

    //Quickly check total particle number
    uint sum=0;
    for (bool S : state)
        if (S)
            ++sum;

    if (sum !=N_particles)
        return {0.0,0.0};

    uint particle_sum= state[0] ? 1:0;
    cx_mat Product = sub_matrices[0+(state[0] ? 1 : 0 )][0];
    for (uint j = 1; j<L; ++j)
    {

        if (state[j])
            ++particle_sum;


        uint INDEXDUM =particle_sum-firstsubcol_totalparticles[j]-top_submatrix_location[j*d+(state[j] ? 1: 0)].second;



        Product = Product*sub_matrices[j*d+(state[j] ? 1 : 0 )][INDEXDUM];
    }

    return Product(0,0);
}

//Stepping left actually makes the matrix at the canonization center right canonical
bool mps::step_left_QR()
{
    if (d!=2)
        throw std::runtime_error("This implementation of step left QR need d=2");

    //cout<<"Stepping left from "<<can_center<<endl;
    if (can_center == -2)
        can_center =L-1;
    if (can_center!= 0 && can_center !=-1)//Nothing to do if everything is all right already
    {


        int j = can_center;

        for (uint k = 0; k<n_subrows[j];++k)
        {
            //Check which matrices have stuff in them, which we can merge

            //INTENTIONAL UNDERFLOW TO DETECT OUT OF RANGE
            uint sig0_id = k-top_submatrix_location[j*d+0].first;
            uint sig1_id = k-top_submatrix_location[j*d+1].first;

            bool sig0_included = sig0_id < sub_matrices[d*j+0].size();
            bool sig1_included = sig1_id < sub_matrices[d*j+1].size();


            uint old_height = subdimension_height[j][k];
            cx_mat L;
            if (sig0_included && sig1_included)
            {
                cx_mat chi=join_rows(
                    sub_matrices[j*d+0][sig0_id],
                    sub_matrices[j*d+1][sig1_id]
                );

                cx_mat Q_dagger;
                cx_mat chi_dagger = chi.t();
                qr_econ(Q_dagger,L,chi_dagger);

                inplace_trans(L);

                inplace_trans(Q_dagger);

                uint width0 = subdimension_height[j+1][top_submatrix_location[j*d+0].second+sig0_id];

                sub_matrices[j*d+0][sig0_id]=Q_dagger.head_cols(width0);
                uint width1 = subdimension_height[j+1][top_submatrix_location[j*d+1].second+sig1_id];

                sub_matrices[j*d+1][sig1_id]=Q_dagger.tail_cols(width1);

                subdimension_height[j][k]=std::min(subdimension_height[j][k],width0+width1);

            }
            else if (sig0_included)
            {
                cx_mat chi_dagger = sub_matrices[j*d+0][sig0_id].t();


                qr_econ(sub_matrices[j*d+0][sig0_id],L,chi_dagger);

                inplace_trans(L);
                inplace_trans(sub_matrices[j*d+0][sig0_id]);

                subdimension_height[j][k]=std::min(subdimension_height[j][k],subdimension_height[j+1][top_submatrix_location[j*d+0].second+sig0_id]);


            }
            else if (sig1_included)
            {
                cx_mat chi_dagger = sub_matrices[j*d+1][sig1_id].t();


                qr_econ(sub_matrices[j*d+1][sig1_id],L,chi_dagger);

                inplace_trans(L);
                inplace_trans(sub_matrices[j*d+1][sig1_id]);


                subdimension_height[j][k]=std::min(subdimension_height[j][k],subdimension_height[j+1][top_submatrix_location[j*d+1].second+sig1_id]);
            }
            else//This should NOT be an option at all
            {
                throw std::runtime_error("QR factorization error, no submatrices at site "+to_string(j)+" row "+to_string(k));
            }
            D_list[j]+=subdimension_height[j][k]-old_height;

            for (uint K = k+1; K<n_subrows[j];++K)
            {
                row_start[j][K]+=subdimension_height[j][k]-old_height;
            }




            uint sig0_prev_id = k-top_submatrix_location[(j-1)*d+0].second;
            uint sig1_prev_id = k-top_submatrix_location[(j-1)*d+1].second;

            if (sig0_prev_id <sub_matrices[d*(j-1)+0].size())
            {
                sub_matrices[(j-1)*d+0][sig0_prev_id]=sub_matrices[(j-1)*d+0][sig0_prev_id]*L;
            }
            if (sig1_prev_id <sub_matrices[d*(j-1)+1].size())
            {
                sub_matrices[(j-1)*d+1][sig1_prev_id]=sub_matrices[(j-1)*d+1][sig1_prev_id]*L;

            }


        }
        can_center = j-1;
    }

    return can_center >0;

}


bool mps::step_right_QR()
{
    if (d!=2)
        throw std::runtime_error("This implementation of step left QR need d=2");

    //cout<<"Stepping left from "<<can_center<<endl;
    if (can_center == -2)
        can_center =0;
    if (can_center!= L-1 && can_center !=L)//Nothing to do if everything is all right already
    {


        int j = can_center;

        for (uint k = 0; k<n_subrows[j+1];++k)
        {
            //Check which matrices have stuff in them, which we can merge

            //INTENTIONAL UNDERFLOW TO DETECT OUT OF RANGE
            uint sig0_id = k-top_submatrix_location[j*d+0].second;
            uint sig1_id = k-top_submatrix_location[j*d+1].second;

            bool sig0_included = sig0_id < sub_matrices[d*j+0].size();
            bool sig1_included = sig1_id < sub_matrices[d*j+1].size();

            uint old_width = subdimension_height[j+1][k];

            cx_mat R;
            if (sig0_included && sig1_included)
            {
                cx_mat psi =join_cols(
                    sub_matrices[j*d+0][sig0_id],
                    sub_matrices[j*d+1][sig1_id]
                );

                cx_mat Q;
                qr_econ(Q,R,psi);

                uint height0 = subdimension_height[j][top_submatrix_location[j*d+0].first+sig0_id];

                sub_matrices[j*d+0][sig0_id]=Q.head_rows(height0);
                uint height1 = subdimension_height[j][top_submatrix_location[j*d+1].first+sig1_id];

                sub_matrices[j*d+1][sig1_id]=Q.tail_rows(height1);

                subdimension_height[j+1][k]=std::min(subdimension_height[j+1][k],height0+height1);
            }
            else if (sig0_included)
            {
                cx_mat psi = sub_matrices[j*d+0][sig0_id];

                qr_econ(sub_matrices[j*d+0][sig0_id],R,psi);

                subdimension_height[j+1][k]=std::min(subdimension_height[j+1][k],subdimension_height[j][top_submatrix_location[j*d+0].first+sig0_id]);


            }
            else if (sig1_included)
            {
                cx_mat psi = sub_matrices[j*d+1][sig1_id];

                qr_econ(sub_matrices[j*d+1][sig1_id],R,psi);

                subdimension_height[j+1][k]=std::min(subdimension_height[j+1][k],subdimension_height[j][top_submatrix_location[j*d+1].first+sig1_id]);
            }
            else//This should NOT be an option at all
            {
                throw std::runtime_error("QR factorization error, no submatrices at site "+to_string(j)+" row "+to_string(k));
            }

            D_list[j+1]+=subdimension_height[j+1][k]-old_width;

            for (uint K = k+1; K<n_subrows[j+1];++K)
            {
                row_start[j+1][K]+=subdimension_height[j+1][k]-old_width;
            }


            uint sig0_next_id = k-top_submatrix_location[(j+1)*d+0].first;
            uint sig1_next_id = k-top_submatrix_location[(j+1)*d+1].first;

            if (sig0_next_id <sub_matrices[d*(j+1)+0].size())
            {
                sub_matrices[(j+1)*d+0][sig0_next_id]=R*sub_matrices[(j+1)*d+0][sig0_next_id];
            }
            if (sig1_next_id <sub_matrices[d*(j+1)+1].size())
            {
                sub_matrices[(j+1)*d+1][sig1_next_id]=R*sub_matrices[(j+1)*d+1][sig1_next_id];
            }
        }
        can_center = j+1;
    }

    return can_center+1 <(int)L;//gcc, please stop these unsigned integer comparison warnings, have an (int) to shut up

}

complex<double> mps::get_coefficient(vector<uint> state) const
{

    //Quickly check total particle number
    uint sum=0;
    for (uint S : state)
        sum+=S;

    if (sum !=N_particles)
        return {0.0,2.0};


    uint particle_sum= state[0];
    cx_mat Product = sub_matrices[0+state[0]][0];
    for (uint j = 1; j<L; ++j)
    {

        particle_sum+= state[j];
        uint INDEXDUM =particle_sum-firstsubcol_totalparticles[j]-top_submatrix_location[j*d+state[j]].second;

        Product = Product*sub_matrices[j*d+(state[j] ? 1 : 0 )][INDEXDUM];
    }

    return Product(0,0);
}

void mps::print_canonization() const
{
    if (can_center == -2)
    {
        cout<<"No canonization has been set"<<endl;
    }
    else
    {
        for (uint j = 0; j<L; ++j)
        {
            cout<<"Site "<<j<<" Compared to "<<can_center<<flush;
            if (j==can_center)
                cout<<" is center"<<endl;
            else if (j<can_center)
            {
                bool approx_ID=true;

                for (uint k = 0; k< n_subrows[j+1]; ++k)
                {


                    bool first=true;
                    cx_mat sum;
                    for (uint sig = 0 ; sig < d; ++sig)
                    {
                        uint this_id = k-top_submatrix_location[j*d+sig].second;



                        if (this_id<sub_matrices[d*j+sig].size())
                        {
                            if (first)
                                    sum=sub_matrices[d*j+sig][this_id].t()*sub_matrices[d*j+sig][this_id];
                            else
                                sum=sum+sub_matrices[d*j+sig][this_id].t()*sub_matrices[d*j+sig][this_id];
                            first =false;
                        }
                    }


                    for (uint n =0; n<subdimension_height[j+1][k]; ++n)
                        for (uint m =0; m<subdimension_height[j+1][k]; ++m)
                        {
                            if (!approx(sum(n,m),n==m ? 1.0 : 0.0))
                                approx_ID=false;
                        }


                }
                if (approx_ID)
                    cout<<" is left"<<endl;
                else
                    cout<<" NOT left"<<endl;
            }
            else if (j>can_center)
            {

                bool approx_ID=true;

                for (uint k = 0; k< n_subrows[j]; ++k)
                {


                    bool first=true;
                    cx_mat sum;
                    for (uint sig = 0 ; sig < d; ++sig)
                    {
                        uint this_id = k-top_submatrix_location[j*d+sig].first;

                        if (this_id<sub_matrices[d*j+sig].size())
                        {
                            if (first)
                                    sum=sub_matrices[d*j+sig][this_id]*sub_matrices[d*j+sig][this_id].t();
                            else
                                sum=sum+sub_matrices[d*j+sig][this_id]*sub_matrices[d*j+sig][this_id].t();
                            first =false;
                        }
                    }


                    for (uint n =0; n<subdimension_height[j][k]; ++n)
                        for (uint m =0; m<subdimension_height[j][k]; ++m)
                        {
                            if (!approx(sum(n,m),n==m ? 1.0 : 0.0))
                                approx_ID=false;
                        }


                }
                if (approx_ID)
                    cout<<" is right"<<endl;
                else
                    cout<<" NOT right"<<endl;
            }
        }
    }
}


//No checks at all, to speed up
complex<double> mps::get(uint sig,uint j, uint k,uint a,uint ap1) const
{
    return sub_matrices[j*d+sig][k](a,ap1);
}

void mps::set(uint sig,uint j, uint k,uint a,uint ap1,complex<double> data)
{
    sub_matrices[j*d+sig][k](a,ap1)=data;
}




//At this j,sig get submatrix k (return false if there is none)
bool mps::sub_matrix_range(uint j, uint sig, uint k, uint& a0, uint& a1, uint& a0p1, uint& a1p1) const
{

   uint row = top_submatrix_location[j*d+sig].first+k;
   uint col = top_submatrix_location[j*d+sig].second+k;

    if (row>=n_subrows[j] || col>=n_subrows[j+1])
        return false;

    a0=row_start[j][row];

    if (row+1==n_subrows[j])
        a1=D_list[j];
    else
        a1=row_start[j][row+1];

    a0p1=row_start[j+1][col];

    if (col+1==n_subrows[j+1])
        a1p1=D_list[j+1];
    else
        a1p1=row_start[j+1][col+1];

    return true;
}

complex<double> mps::overlap(const mps& that) const
{
    //Only works if we have the same block-structure!, then, in that case this is extremely efficient, the middle matrix in the overlap simply becomes a block-diagonal matrix

    if (n_subrows[1]!=that.n_subrows[1])
        throw std::runtime_error("Block structure does not match!");

    vector<cx_mat> center_blocks(n_subrows[1]);

    for (uint k = 0; k<n_subrows[1]; ++k)
    {
        center_blocks[k]=cx_mat(subdimension_height[1][k],that.subdimension_height[1][k]);


    }

    for (uint sig = 0; sig<d; ++sig)
    {
        for (uint i = 0; i<sub_matrices[0*d+sig].size(); ++i)
        {

            //This diagonal element is in this row/column in the combined matrix
            uint k =top_submatrix_location[0*d+sig].second+i;
            center_blocks[k]+=sub_matrices[0*d+sig][i].t()*that.sub_matrices[0*d+sig][i];

        }
    }

    //Now loop through all sites:

    for (uint j=1; j<L; ++j)
    {
        vector<cx_mat> new_center_blocks(n_subrows[1+j]);
        for (uint k = 0; k<n_subrows[j+1]; ++k)
        {
            new_center_blocks[k]=cx_mat(subdimension_height[j+1][k],that.subdimension_height[j+1][k]);
        }

        for (uint sig = 0; sig<d; ++sig)
        {
            for (uint i = 0; i<sub_matrices[j*d+sig].size(); ++i)
            {
                //This diagonal element is in this row/column in the combined matrix
                uint row =top_submatrix_location[j*d+sig].second+i;
                uint col=top_submatrix_location[j*d+sig].first+i;



                new_center_blocks[row]+=((sub_matrices[j*d+sig][i].t())*center_blocks[col])*that.sub_matrices[j*d+sig][i];
            }
        }

        center_blocks = std::move(new_center_blocks);
    }

    return center_blocks[0](0,0);
}


//We need the largest integer we can get when working with factorials
uint64_t factorial(uint64_t n)
{
    uint64_t out = 1;
    for (;n>0;--n)
        out*=n;
    return out;
}

//There are likely smarter ways of doing this, but this method does at least not try to store all of n! which is likely to overflow. We do still need to calculate and store whichever is smaller of m! or (n-m)! which is up to (n/2)! overflow is still likely so I will just use uint64_t
uint64_t bin_coef(uint64_t n, uint64_t m)
{
    if (n>=m && m>=0)
    {
        //Overflow is extremely likely if we use our formula
        //factorial(n)/(factorial(m)*factorial(n-m));
        //But in practice we are taking (n)*(n-1)*(n-2)*... (m+1) or (n-m+1), whichever is larger, divided by (n-m)! or m! whichever is smaller, lets label the larger m

        uint64_t M;//Whichever is larger of m and n-m; (this is really just bin_coef(n,m)=bin_coef(n,n-m) )

        if (2*m > n)
        {
            M=m;
            m=n-m;//Use m to store whichever is smaller
        }
        else
            M=n-m;

        uint64_t out = 1;
        for (uint64_t i = n; i>M; --i)
            out*=i;//Calculate (n)*(n-1)*(n-2)*... (m+1) or (n-m+1), whichever is larger
        out/=factorial(m);//Divide with whatever is smaller

        return out;
    }
    else
        return 0;
}

