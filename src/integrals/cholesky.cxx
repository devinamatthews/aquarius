#include "cholesky.hpp"

using namespace aquarius::tensor;
using namespace aquarius::input;
using namespace aquarius::task;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace integrals
{

template <typename T>
CholeskyIntegrals<T>::CholeskyIntegrals(const Arena& arena, const Context& ctx, Config& config, const Molecule& molecule)
: Distributed(arena),
  molecule(molecule),
  ctx(ctx),
  nvec(0),
  shells(molecule.getShellsBegin(),molecule.getShellsEnd()),
  delta(config.get<T>("delta")),
  cond(config.get<T>("cond_max"))
{
    nfunc = 0;
    for (int i = 0;i < shells.size();i++) nfunc += shells[i].getNFunc()*shells[i].getNContr();
    ndiag = nfunc * (nfunc + 1) / 2;
    nblock = shells.size() * (shells.size() + 1) / 2;
    decompose();
}

template <typename T>
void CholeskyIntegrals<T>::test()
{
    const PointGroup& group = molecule.getGroup();
    SymmetryBlockedTensor<T> LD("LD", this->arena, group, 3, {{nfunc},{nfunc},{ndiag}}, {SY,NS,NS}, false);
    SymmetryBlockedTensor<T> ints("V", this->arena, group, 4, {{nfunc},{nfunc},{nfunc},{nfunc}}, {NS,NS,NS,NS}, false);

    LD["pqJ"] = (*L)["pqJ"]*(*D)["J"];
    ints["pqrs"] = (*L)["pqJ"]*LD["rsJ"];

    //cout << "ints:" << endl;
    //ints.print(stdout);
    //cout << endl;

    T err = 0;
    for (int i = 0, m = 0;i < shells.size();i++)
    {
        for (int j = 0, n = 0;j <= i;j++)
        {
            for (int k = 0, o = 0;k < shells.size();k++)
            {
                for (int l = 0, p = 0;l <= k;l++)
                {
                    int ni = shells[i].getNFunc()*shells[i].getNContr();
                    int nj = shells[j].getNFunc()*shells[j].getNContr();
                    int nk = shells[k].getNFunc()*shells[k].getNContr();
                    int nl = shells[l].getNFunc()*shells[l].getNContr();

                    DenseTensor<T> local_ints("tmp", 4, {ni,nj,nk,nl});

                    vector<tkv_pair<T>> pairs(ni*nj*nk*nl);

                    for (int a = 0, q = 0;a < ni;a++)
                    {
                        for (int b = 0;b < nj;b++)
                        {
                            for (int c = 0;c < nk;c++)
                            {
                                for (int d = 0;d < nl;d++)
                                {
                                    pairs[q++].k = (((d+p)*nfunc+(c+o))*nfunc+(b+n))*nfunc+(a+m);
                                }
                            }
                        }
                    }

                    if (arena.rank == 0)
                    {
                        ints(0).getRemoteData(pairs);

                        T* local_data = local_ints.getData();
                        for (int q = 0;q < ni*nj*nk*nl;q++)
                        {
                            int d = pairs[q].k/(nfunc*nfunc*nfunc) - p;
                            pairs[q].k -= (d+p)*nfunc*nfunc*nfunc;
                            int c = pairs[q].k/(nfunc*nfunc) - o;
                            pairs[q].k -= (c+o)*nfunc*nfunc;
                            int b = pairs[q].k/(nfunc) - n;
                            pairs[q].k -= (b+n)*nfunc;
                            int a = pairs[q].k - m;

                            local_data[((d*nk+c)*nj+b)*ni+a] = pairs[q].d;
                        }

                        //printf("%d %d %d %d\n", i, j, k, l);
                        err += testBlock(local_ints, shells[i], shells[j], shells[k], shells[l]);
                    }
                    else
                    {
                        ints(0).getRemoteData();
                    }

                    p += shells[l].getNFunc()*shells[l].getNContr();
                }
                o += shells[k].getNFunc()*shells[k].getNContr();
            }
            n += shells[j].getNFunc()*shells[j].getNContr();
        }
        m += shells[i].getNFunc()*shells[i].getNContr();
    }

    if (arena.rank == 0)
    {
        err = sqrt(err / ndiag / ndiag);
        printf("RMS error: %.15e\n\n", err);
    }
}

template <typename T>
void CholeskyIntegrals<T>::decompose()
{
    diag_elem_t* diag = new diag_elem_t[ndiag];
    for (int j = 0, elem = 0;j < shells.size();j++)
    {
        for (int i = 0;i <= j;i++)
        {
            elem += getDiagonalBlock(shells[i], shells[j], diag+elem);
        }
    }

    int* block_start = new int[nblock];
    int* block_size = new int[nblock];
    sortBlocks(diag, block_start, block_size);

    int nblock_local = 0;
    for (int elem = 0;;)
    {
        int block = arena.rank+nblock_local*arena.size;

        if (block >= nblock) break;

        copy(diag+block_start[block], diag+block_start[block]+block_size[block], diag+elem);
        block_start[nblock_local] = elem;
        block_size[nblock_local] = block_size[block];
        elem += block_size[nblock_local];

        nblock_local++;
    }

    T** block_data = new T*[nblock_local];
    int max_block_size = 0;
    for (int block = 0;block < nblock_local;block++)
    {
        max_block_size = max(max_block_size, block_size[block]);
        block_data[block] = new T[block_size[block]*ndiag];
        fill(block_data[block], block_data[block]+block_size[block]*ndiag, 0.0);
    }

    arena.Allreduce(&max_block_size, 1, MPI_MAX);

    //for (l = 0;l < ndiag;l++)
    //{
    //    if (diag[l].elem < delta * delta / crit) diag[l].elem = 0;
    //}

    T* D = new T[ndiag];
    T* tmp_block_data = new T[ndiag*max_block_size];
    diag_elem_t* tmp_diag = new diag_elem_t[max_block_size];
    for (nvec = 0;;)
    {
        int converged = 1;
        T local_max = -2;
        int max_block = 0;
        for (int block = 0;block < nblock_local;block++)
        {
            T max_elem;
            //printf("checking block: %d\n", block);
            converged = isBlockConverged(diag+block_start[block], block_size[block], max_elem) && converged;
            if (max_elem > local_max)
            {
                local_max = max_elem;
                max_block = block;
            }
        }

        //printf("max block: %d\n", max_block);
        //printf("max elem: %e\n", local_max);

        arena.Allreduce(&converged, 1, MPI_BAND);
        if (converged) break;

        vector<T> maxes(arena.size);
        arena.Allgather(&local_max, maxes);

        int pmax = 0;
        T global_max = 0;
        for (int p = 0;p < arena.size;p++)
        {
            if (maxes[p] > global_max)
            {
                global_max = maxes[p];
                pmax = p;
            }
        }

        int old_rank = nvec;
        int nactive;

        if (pmax == arena.rank)
        {
            //cout << "Decomposing block " << rank+max_block*nproc << endl;

            decomposeBlock(block_size[max_block], block_data[max_block], D, diag+block_start[max_block]);

            nactive = collectActiveRows(block_size[max_block], block_data[max_block], diag+block_start[max_block],
                                        tmp_block_data, tmp_diag);
        }

        arena.Bcast(&nvec, 1, pmax);
        if (nvec == old_rank) continue;

        arena.Bcast(&nactive, 1, pmax);
        arena.Bcast(tmp_block_data, nactive*ndiag, pmax);
        arena.Bcast(tmp_diag, nactive*sizeof(diag_elem_t), pmax, MPI_BYTE);
        arena.Bcast(D+old_rank, nvec-old_rank, pmax);

        for (int next_block = 0;next_block < nblock_local;next_block++)
        {
            if (next_block == max_block && arena.rank == pmax) continue;

            updateBlock(old_rank, block_size[next_block], block_data[next_block], diag+block_start[next_block],
                                  nactive,                tmp_block_data,         tmp_diag,                     D);
        }
    }

    assert(nvec > 0);

    printf("rank: full, partial: %d %d\n\n", ndiag, nvec);

    for (int block = 0;block < nblock_local;block++)
        resortBlock(block_size[block], block_data[block], diag+block_start[block], tmp_block_data);

    delete[] tmp_block_data;
    delete[] tmp_diag;

    const PointGroup& group = molecule.getGroup();

    this->D = new SymmetryBlockedTensor<T>("D", this->arena, group, 1, {{nvec}}, {NS}, false);
    this->L = new SymmetryBlockedTensor<T>("L", this->arena, group, 3, {{nfunc},{nfunc},{nvec}}, {SY,NS,NS}, false);

    if (arena.rank == 0)
    {
        vector<tkv_pair<T>> pairs(nvec);
        for (int i = 0;i < nvec;i++)
        {
            pairs[i].k = i;
            pairs[i].d = D[i];
        }
        (*this->D)(0).writeRemoteData(pairs);
    }
    else
    {
        (*this->D)(0).writeRemoteData();
    }

    vector<vector<int>> idx = Shell::setupIndices(ctx, molecule);

    vector<tkv_pair<T>> pairs;
    for (int block = 0;block < nblock_local;block++)
    {
        int i = diag[block_start[block]].shelli;
        int j = diag[block_start[block]].shellj;

        //printf("\nblock %d\n", block);
        //print_matrix(block_data[block], block_size[block], rank);

        int elem = 0;
        for (int f = 0;f < shells[j].getNFunc();f++)
        {
            for (int n = 0;n < shells[j].getNContr();n++)
            {
                for (int e = 0;e < shells[i].getNFunc();e++)
                {
                    for (int m = 0;m < shells[i].getNContr();m++)
                    {
                        if (i == j &&
                            e*shells[i].getNContr()+m >
                            f*shells[j].getNContr()+n) break;

                        int o = shells[i].getIndex(ctx, idx[i], e, m, 0);
                        int p = shells[j].getIndex(ctx, idx[j], f, n, 0);
                        //key idx = max(o,p)*(max(o,p)+1)/2 + min(o,p);
                        key idx = max(o,p)*nfunc + min(o,p);
                        //printf("%d %d %d %d\n", i, j, shells[i].getIdx()[0], shells[j].getIdx()[0]);
                        for (int r = 0;r < nvec;r++)
                        {
                            //printf("%d %d %d %d %d %d %d %d %d %d\n", rank, idx, o, p, i, j, m, n, r, block);
                            pairs.push_back(tkv_pair<T>(idx, block_data[block][elem*nvec+r]));
                            //idx += ndiag;
                            idx += nfunc*nfunc;
                        }

                        elem++;
                    }
                }
            }
        }
    }
    (*this->L)(0).writeRemoteData(pairs);

    //cout << "D:" << endl;
    //this->D->print(cout);
    //cout << endl;

    //cout << "L:" << endl;
    //this->L->print(cout);
    //cout << endl;

    for (int i = 0;i < nblock_local;i++) delete[] block_data[i];
    delete[] D;
    delete[] diag;
    delete[] block_data;
    delete[] block_size;
    delete[] block_start;
}

template <typename T>
void CholeskyIntegrals<T>::resortBlock(const int block_size, T* L, diag_elem_t* diag, T* tmp)
{
    for (int elem = 0;elem < block_size;elem++)
    {
        dcopy(nvec, L+elem*ndiag, 1, tmp+diag[elem].idx*nvec, 1);
    }

    dcopy(nvec*block_size, tmp, 1, L, 1);
}

template <typename T>
int CholeskyIntegrals<T>::collectActiveRows(const int block_size, const T* L, diag_elem_t* diag,
                      T* L_active, diag_elem_t* diag_active)
{
    int cur = 0;

    for (int elem = 0;elem < block_size;elem++)
    {
        if (diag[elem].status == ACTIVE)
        {
            //printf("finishing row: %d\n", elem);
            diag[elem].status = DONE;
            diag_active[cur] = diag[elem];
            dcopy(nvec, L+elem*ndiag, 1, L_active+cur*ndiag, 1);
            cur++;
        }
    }

    return cur;
}

template <typename T>
bool CholeskyIntegrals<T>::isBlockConverged(const diag_elem_t* diag, int block_size, T& max_elem)
{
    bool converged = true;
    max_elem = -1;

    for (int elem = 0;elem < block_size;elem++)
    {
        if (fabs(diag[elem].elem) > delta && diag[elem].status == TODO)
        {
            max_elem = max(max_elem, fabs(diag[elem].elem));
            //printf("row %d isn't done yet: %e %e\n", elem, max_elem, fabs(diag[elem].elem));
            converged &= false;
        }
    }

    return converged;
}

template <typename T>
void CholeskyIntegrals<T>::sortBlocks(diag_elem_t* diag, int* block_start, int* block_size)
{
    sort(diag, diag+ndiag);

    int block = 0;
    for (int cur = 0;cur < ndiag;cur++)
    {
        assert(block < nblock);

        block_start[block] = cur;

        for (int end = cur+1;end < ndiag;end++)
        {
            if (diag[cur].shelli == diag[end].shelli &&
                diag[cur].shellj == diag[end].shellj)
            {
                rotate(diag+(++cur), diag+end, diag+end+1);
            }
        }

        block_size[block] = (cur+1)-block_start[block];
        block++;
    }

    assert(block == nblock);
}

template <typename T>
void CholeskyIntegrals<T>::getShellOffsets(
    const Shell& a, const Shell& b, const Shell& c, const Shell& d,
    size_t& controffa, size_t& funcoffa, size_t& controffb, size_t& funcoffb,
    size_t& controffc, size_t& funcoffc, size_t& controffd, size_t& funcoffd)
{
    controffa = 0;
    controffb = 0;
    controffc = 0;
    controffd = 0;
    funcoffa = 0;
    funcoffb = 0;
    funcoffc = 0;
    funcoffd = 0;

    size_t ncab = a.getNContr()*b.getNContr();
    size_t nccd = c.getNContr()*d.getNContr();
    size_t nfab = a.getNFunc()*b.getNFunc();
    size_t nfcd = c.getNFunc()*d.getNFunc();

    /*
    if (&context.getA() == &a || &context.getA() == &b)
    {
        if (&context.getA() == &a)
        {
    */
            controffa += 1;
            controffb += a.getNContr();
            funcoffa += ncab*nccd;
            funcoffb += ncab*nccd*a.getNFunc();
   /*
        }
        else
        {
            controffa += b.getNContr();
            controffb += 1;
            funcoffa += ncab*nccd*b.getNFunc();
            funcoffb += ncab*nccd;
        }

        if (&context.getC() == &c)
        {
    */
            controffc += ncab;
            controffd += ncab*c.getNContr();
            funcoffc += ncab*nccd*nfab;
            funcoffd += ncab*nccd*nfab*c.getNFunc();
    /*
        }
        else
        {
            controffc += ncab*d.getNContr();
            controffd += ncab;
            funcoffc += ncab*nccd*nfab*d.getNFunc();
            funcoffd += ncab*nccd*nfab;
        }
    }
    else
    {
        if (&context.getA() == &c)
        {
            controffc += 1;
            controffd += c.getNContr();
            funcoffc += ncab*nccd;
            funcoffd += ncab*nccd*c.getNFunc();
        }
        else
        {
            controffc += d.getNContr();
            controffd += 1;
            funcoffc += ncab*nccd*d.getNFunc();
            funcoffd += ncab*nccd;
        }

        if (&context.getC() == &a)
        {
            controffa += nccd;
            controffb += nccd*a.getNContr();
            funcoffa += ncab*nccd*nfcd;
            funcoffb += ncab*nccd*nfcd*a.getNFunc();
        }
        else
        {
            controffa += nccd*b.getNContr();
            controffb += nccd;
            funcoffa += ncab*nccd*nfcd*b.getNFunc();
            funcoffb += ncab*nccd*nfcd;
        }
    }
    */
}

template <typename T>
int CholeskyIntegrals<T>::getDiagonalBlock(const Shell& a, const Shell& b, diag_elem_t* diag)
{
    TwoElectronIntegrals eri(a, b, a, b, ERIEvaluator());
    const T* intbuf = eri.getIntegrals();

    size_t controffi;
    size_t funcoffi;
    size_t controffj;
    size_t funcoffj;

    getShellOffsets(a, b, a, b,
                    controffi, funcoffi, controffj, funcoffj,
                    controffi, funcoffi, controffj, funcoffj);

    int elem = 0;
    for (int n = 0;n < b.getNFunc();n++)
    {
        for (int p = 0;p < b.getNContr();p++)
        {
            for (int m = 0;m < a.getNFunc();m++)
            {
                for (int o = 0;o < a.getNContr();o++)
                {
                    if (&a == &b && (m*a.getNContr()+o) > (n*b.getNContr()+p)) continue;

                    diag[elem].funci = m;
                    diag[elem].contri = o;
                    diag[elem].funcj = n;
                    diag[elem].contrj = p;
                    diag[elem].elem = intbuf[m*funcoffi+o*controffi+
                                             n*funcoffj+p*controffj];
                    diag[elem].shelli = (int)(&a-&shells[0]);
                    diag[elem].shellj = (int)(&b-&shells[0]);
                    diag[elem].idx = elem;
                    diag[elem].status = TODO;

                    //crit = max(crit, diag[elem].elem);
                    elem++;
                }
            }
        }
    }

    return elem;
}

template <typename T>
void CholeskyIntegrals<T>::decomposeBlock(int block_size, T* L_, T* D, diag_elem_t* diag)
{
    //printf("starting block: %d %f %d\n", block_size, diag[0].elem, rank);

    T (*L)[ndiag] = (T(*)[ndiag])L_;

    bool found = false;
    T max_elem = 0;
    for (int elem = 0;elem < block_size;elem++)
    {
        if (diag[elem].status != TODO) continue;
        max_elem = max(max_elem, fabs(diag[elem].elem));
        if (fabs(diag[elem].elem) > delta) found = true;
    }

    T crit = max(delta, max_elem/cond);

    if (!found) return;

    TwoElectronIntegrals eri(shells[diag[0].shelli], shells[diag[0].shellj],
                             shells[diag[0].shelli], shells[diag[0].shellj], ERIEvaluator());
    const T* intbuf = eri.getIntegrals();

    size_t controffi;
    size_t funcoffi;
    size_t controffj;
    size_t funcoffj;
    size_t controffk;
    size_t funcoffk;
    size_t controffl;
    size_t funcoffl;

    getShellOffsets(shells[diag[0].shelli], shells[diag[0].shellj],
                    shells[diag[0].shelli], shells[diag[0].shellj],
                    controffi, funcoffi, controffj, funcoffj,
                    controffk, funcoffk, controffl, funcoffl);

    for (int elem = 0;elem < block_size;elem++)
    {
        //printf("activating row: %d\n", elem);
        if (fabs(diag[elem].elem) < crit || diag[elem].status == DONE) continue;

        diag[elem].status = ACTIVE;

        //printf("updating L[%d][%d] (out of %d %d, addr=0x%x)\n", rank, elem, ndiag, block_size, L_);
        L[elem][nvec] = 1.0;
        D[nvec] = diag[elem].elem;

        #pragma omp parallel for schedule(dynamic) default(shared)
        for (int row = 0;row < block_size;row++)
        {
            if (diag[row].status != TODO) continue;

            //printf("updating L[%d][%d]\n", rank, row);
            L[row][nvec] = intbuf[diag[elem].contri*controffi+diag[elem].funci*funcoffi+
                                  diag[elem].contrj*controffj+diag[elem].funcj*funcoffj+
                                  diag[ row].contri*controffk+diag[ row].funci*funcoffk+
                                  diag[ row].contrj*controffl+diag[ row].funcj*funcoffl];
            for (int col = 0;col < nvec;col++)
                L[row][nvec] -= D[col] * L[row][col] * L[elem][col];
            L[row][nvec] /= D[nvec];
            diag[row].elem -= D[nvec] * L[row][nvec] * L[row][nvec];
        }

        nvec++;
    }
}

template <typename T>
void CholeskyIntegrals<T>::updateBlock(int old_rank, int block_size_i, T* L_i_, diag_elem_t* diag_i,
                                       int block_size_j, T* L_j_, diag_elem_t* diag_j, T* D)
{
    T (*L_i)[ndiag] = (T(*)[ndiag])L_i_;
    T (*L_j)[ndiag] = (T(*)[ndiag])L_j_;

    bool found = false;
    for (int elem = 0;elem < block_size_i;elem++)
    {
        if (diag_i[elem].status == TODO) found = true;
    }

    if (!found) return;

    //printf("subblock: %d %d\n", l, diag[l].nblock);

    TwoElectronIntegrals eri(shells[diag_j[0].shelli], shells[diag_j[0].shellj],
                             shells[diag_i[0].shelli], shells[diag_i[0].shellj], ERIEvaluator());
    const T* intbuf = eri.getIntegrals();

    size_t controffii;
    size_t funcoffii;
    size_t controffij;
    size_t funcoffij;
    size_t controffji;
    size_t funcoffji;
    size_t controffjj;
    size_t funcoffjj;

    getShellOffsets(shells[diag_j[0].shelli], shells[diag_j[0].shellj],
                    shells[diag_i[0].shelli], shells[diag_i[0].shellj],
                    controffji, funcoffji, controffjj, funcoffjj,
                    controffii, funcoffii, controffij, funcoffij);

    #pragma omp parallel for schedule(dynamic) default(shared)
    for (int elem_i = 0;elem_i < block_size_i;elem_i++)
    {
        if (diag_i[elem_i].status != TODO) continue;

        int cur_rank = old_rank;
        for (int elem_j = 0;elem_j < block_size_j;elem_j++)
        {
            //printf("updating L[%d][%d]\n", cur_rank, elem_i);
            L_i[elem_i][cur_rank] = intbuf[diag_j[elem_j].contri*controffji+diag_j[elem_j].funci*funcoffji+
                                           diag_j[elem_j].contrj*controffjj+diag_j[elem_j].funcj*funcoffjj+
                                           diag_i[elem_i].contri*controffii+diag_i[elem_i].funci*funcoffii+
                                           diag_i[elem_i].contrj*controffij+diag_i[elem_i].funcj*funcoffij];
            for (int col = 0;col < cur_rank;col++)
                L_i[elem_i][cur_rank] -= D[col] * L_i[elem_i][col] * L_j[cur_rank-old_rank][col];
            L_i[elem_i][cur_rank] /= D[cur_rank];
            diag_i[elem_i].elem -= D[cur_rank] * L_i[elem_i][cur_rank] * L_i[elem_i][cur_rank];

            cur_rank++;
        }
    }
}

template <typename T>
T CholeskyIntegrals<T>::testBlock(const DenseTensor<T>& block, const Shell& a, const Shell& b, const Shell& c, const Shell& d)
{
    DenseTensor<T> packed("tmp", 4, block.getLengths(), false);
    const T* ints = packed.getData();

    packed = block;

    TwoElectronIntegrals eri(a, b, c, d, ERIEvaluator());
    const T* intbuf = eri.getIntegrals();

    size_t controffa;
    size_t funcoffa;
    size_t controffb;
    size_t funcoffb;
    size_t controffc;
    size_t funcoffc;
    size_t controffd;
    size_t funcoffd;

    getShellOffsets(a, b, c, d,
                    controffa, funcoffa, controffb, funcoffb,
                    controffc, funcoffc, controffd, funcoffd);

    double err = 0;
    int q = 0;
    for (int l = 0;l < d.getNFunc();l++)
    {
        for (int p = 0;p < d.getNContr();p++)
        {
            for (int k = 0;k < c.getNFunc();k++)
            {
                for (int o = 0;o < c.getNContr();o++)
                {
                    for (int j = 0;j < b.getNFunc();j++)
                    {
                        for (int n = 0;n < b.getNContr();n++)
                        {
                            for (int i = 0;i < a.getNFunc();i++)
                            {
                                for (int m = 0;m < a.getNContr();m++)
                                {
                                    //if ((&a == &b && (j > i || (j == i && n > m))) ||
                                    //    (&c == &d && (l > k || (l == k && p > o))))
                                    //{
                                    //    q++;
                                    //    continue;
                                    //}

                                    err += pow(ints[q] - intbuf[m*controffa+i*funcoffa+
                                                                n*controffb+j*funcoffb+
                                                                o*controffc+k*funcoffc+
                                                                p*controffd+l*funcoffd], 2);

                                    if (fabs(ints[q] - intbuf[m*controffa+i*funcoffa+
                                                              n*controffb+j*funcoffb+
                                                              o*controffc+k*funcoffc+
                                                              p*controffd+l*funcoffd]) > 1e-10)
                                    {
                                    ///*
                                    printf("%d %d %d %d %d %d %d %d: %.12f %.12f\n",
                                           i, j, k, l, m, n, o, p, ints[q],
                                           intbuf[m*controffa+i*funcoffa+
                                                  n*controffb+j*funcoffb+
                                                  o*controffc+k*funcoffc+
                                                  p*controffd+l*funcoffd]);
                                    //*/
                                    }

                                    q++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return err;
}

INSTANTIATE_SPECIALIZATIONS(CholeskyIntegrals);

}
}
