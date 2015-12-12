#include "fcidump.hpp"

using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::symmetry;

using namespace std::regex_constants;

static constexpr int BUFFER_SIZE = 10000000;

namespace aquarius
{
namespace op
{

template <typename T>
void writeIntegrals(vector<kv_pair>& buf, CTFTensor<T>& tensor)
{
    cout << "dump " << tensor.getLengths() << endl;
    tensor.writeRemoteData(buf);
    buf.clear();
}

template <typename T>
void writeIntegral(int64_t key, double val, vector<kv_pair>& buf, CTFTensor<T>& tensor)
{
    //cout << key << " " << val << endl;
    buf.emplace_back(key, val);
    if (buf.size() == BUFFER_SIZE) writeIntegrals(buf, tensor);
}

template <typename T>
FCIDUMP<T>::FCIDUMP(const string& name, Config& config)
: Task(name, config), path(config.get<string>("filename"))
{
    addProduct("moints", "H");
}

template <typename T>
bool FCIDUMP<T>::run(TaskDAG& dag, const Arena& arena)
{
    ifstream ifs(path);
    string header, line;
    smatch m;
    int norb, nelec, no, nv;

    while (getline(ifs, line))
    {
        header = header + " " + line;
        if (regex_search(line, regex("(/|[$&]END)", icase))) break;
    }

    regex_search(header, m, regex("NORB\\s*=\\s*([0-9]+)", icase));
    istringstream(m[1]) >> norb;

    regex_search(header, m, regex("NELEC\\s*=\\s*([0-9]+)", icase));
    istringstream(m[1]) >> nelec;

    no = nelec/2;
    nv = norb-no;

    this->log(arena) << "There are " << no << " occupied and " << nv << " virtual orbitals." << endl;

    Space occ(PointGroup::C1(), {no}, {no});
    Space vrt(PointGroup::C1(), {nv}, {nv});

    auto& H = this->put("H", new TwoElectronOperator<T>("H", arena, occ, vrt));

    double escf = 0;

    matrix<double> fij(no, no);
    matrix<double> fia(no, nv);
    matrix<double> fai(nv, no);
    matrix<double> fab(nv, nv);
    vector<kv_pair> ijkl_buf; ijkl_buf.reserve(BUFFER_SIZE);
    vector<kv_pair> aijk_buf; aijk_buf.reserve(BUFFER_SIZE);
    vector<kv_pair> abij_buf; abij_buf.reserve(BUFFER_SIZE);
    vector<kv_pair> aibj_buf; aibj_buf.reserve(BUFFER_SIZE);
    vector<kv_pair> abci_buf; abci_buf.reserve(BUFFER_SIZE);
    vector<kv_pair> abcd_buf; abcd_buf.reserve(BUFFER_SIZE);

    CTFTensor<T>& VIJKL = H.getIJKL()({0,1},{0,1})({0,0,0,0});
    CTFTensor<T>& VAIJK = H.getAIJK()({1,0},{0,1})({0,0,0,0});
    CTFTensor<T>& VABIJ = H.getABIJ()({1,0},{0,1})({0,0,0,0});
    CTFTensor<T>& VAIBJ = H.getAIBJ()({1,0},{1,0})({0,0,0,0});
    CTFTensor<T>& VABCI = H.getABCI()({1,0},{1,0})({0,0,0,0});
    CTFTensor<T>& VABCD = H.getABCD()({1,0},{1,0})({0,0,0,0});

    //vector<int> orbmap{0, 1, 2, 5, 6, 7, 8, 11, 3, 9, 4, 10};

    while (getline(ifs, line))
    {
        double val;
        int64_t p, q, r, s;

        istringstream(line) >> val >> p >> q >> r >> s;
        //p = orbmap[p];
        //q = orbmap[q];
        //r = orbmap[r];
        //s = orbmap[s];

        if (p == 0)
        {
            escf += val;
        }
        else if (q == 0)
        {
            continue;
        }
        else if (r == 0)
        {
            bool p_is_vrt = --p >= no; if (p_is_vrt) p -= no;
            bool q_is_vrt = --q >= no; if (q_is_vrt) q -= no;

            if (p_is_vrt)
            {
                if (q_is_vrt)
                {
                    fab[p][q] += val;
                    if (p != q) fab[q][p] += val;
                }
                else
                {
                    fai[p][q] += val;
                    fia[q][p] += val;
                }
            }
            else
            {
                if (q_is_vrt)
                {
                    //hermitian
                    //fia[p][q] += val;
                }
                else
                {
                    if (p == q) escf += 2*val;
                    fij[p][q] += val;
                    if (p != q) fij[q][p] += val;
                }
            }
        }
        else
        {
            val /= 2;

            /*
             * Switch to <pq|rs> with p>=r, q>=s, pr>=qs.
             */
            if (p < q) swap(p, q);
            if (r < s) swap(r, s);
            if (p < r || (p == r && q < s))
            {
                continue;
                swap(p, r);
                swap(q, s);
            }
            swap(q, r);

            bool p_is_vrt = --p >= no; if (p_is_vrt) p -= no;
            bool q_is_vrt = --q >= no; if (q_is_vrt) q -= no;
            bool r_is_vrt = --r >= no; if (r_is_vrt) r -= no;
            bool s_is_vrt = --s >= no; if (s_is_vrt) s -= no;

            bool pr_eq_qs = min(p,r) == min(q,s) && max(p,r) == max(q,s);

            for (int pr = 0;pr < 2;pr++)
            {
                for (int qs = 0;qs < 2;qs++)
                {
                    for (int prqs = 0;prqs < 2;prqs++)
                    {
                        if (r_is_vrt)
                        {
                            if (s_is_vrt)
                            {
                                /*
                                 * VVVV
                                 */
                                writeIntegral(((s*nv+r)*nv+q)*nv+p, val, abcd_buf, VABCD);
                            }
                            else if (q_is_vrt)
                            {
                                /*
                                 * VVVO
                                 */
                                writeIntegral(((s*nv+r)*nv+q)*nv+p, val, abci_buf, VABCI);
                            }
                            else
                            {
                                /*
                                 * VOVO
                                 */
                                if (q == s) fab[p][r] += 2*val;
                                writeIntegral(((s*nv+r)*no+q)*nv+p, val, aibj_buf, VAIBJ);
                            }
                        }
                        else if (p_is_vrt)
                        {
                            if (s_is_vrt)
                            {
                                /*
                                 * VVOV
                                 */
                                writeIntegral(((r*nv+s)*nv+p)*nv+q, val, abci_buf, VABCI);
                            }
                            else if (q_is_vrt)
                            {
                                /*
                                 * VVOO
                                 */
                                if (r == s) fab[p][q] -= val;
                                writeIntegral(((s*no+r)*nv+q)*nv+p, val, abij_buf, VABIJ);
                            }
                            else
                            {
                                /*
                                 * VOOO
                                 */
                                if (q == s) fai[p][r] += 2*val;
                                if (q == s) fia[r][p] += 2*val;
                                if (q == r) fai[p][s] -= val;
                                if (q == r) fia[s][p] -= val;
                                writeIntegral(((s*no+r)*no+q)*nv+p, val, aijk_buf, VAIJK);
                            }
                        }
                        else
                        {
                            if (s_is_vrt)
                            {
                                /*
                                 * OVOV
                                 */
                                abort();
                            }
                            else if (q_is_vrt)
                            {
                                /*
                                 * OVOO
                                 */
                                abort();
                            }
                            else
                            {
                                /*
                                 * OOOO
                                 */
                                if (q == s) fij[p][r] += 2*val;
                                if (q == r) fij[p][s] -= val;
                                if (p == r && q == s) escf += 2*val;
                                if (p == s && q == r) escf -= val;
                                writeIntegral(((s*no+r)*no+q)*no+p, val, ijkl_buf, VIJKL);
                            }
                        }

                        if (pr_eq_qs || p_is_vrt != q_is_vrt || r_is_vrt != s_is_vrt) break;
                        swap(p, q);
                        swap(r, s);
                    }
                    if (q == s || q_is_vrt != s_is_vrt) break;
                    swap(q, s);
                }
                if (p == r || p_is_vrt != r_is_vrt) break;
                swap(p, r);
            }
        }
    }

    double norm_ij = 0;
    vector<kv_pair> ij_buf;
    for (int i = 0;i < no;i++)
    {
        for (int j = 0;j < no;j++)
        {
            if (i != j) norm_ij = max(norm_ij, fabs(fij[i][j]));
            ij_buf.emplace_back(i+j*no, fij[i][j]);
        }
    }

    double norm_ia = 0;
    vector<kv_pair> ia_buf;
    for (int i = 0;i < no;i++)
    {
        for (int a = 0;a < nv;a++)
        {
            norm_ia = max(norm_ia, fabs(fia[i][a]));
            ia_buf.emplace_back(i+a*no, fia[i][a]);
        }
    }

    double norm_ai = 0;
    vector<kv_pair> ai_buf;
    for (int a = 0;a < nv;a++)
    {
        for (int i = 0;i < no;i++)
        {
            norm_ai = max(norm_ai, fabs(fai[a][i]));
            ai_buf.emplace_back(a+i*nv, fai[a][i]);
        }
    }

    double norm_ab = 0;
    vector<kv_pair> ab_buf;
    for (int a = 0;a < nv;a++)
    {
        for (int b = 0;b < nv;b++)
        {
            if (a != b) norm_ab = max(norm_ab, fabs(fab[a][b]));
            ab_buf.emplace_back(a+b*nv, fab[a][b]);
        }
    }

    cout << "E(SCF): " << printToAccuracy(escf, 1e-12) << endl;
    cout << "norm IJ: " << norm_ij << endl;
    cout << "norm IA: " << norm_ia << endl;
    cout << "norm AI: " << norm_ai << endl;
    cout << "norm AB: " << norm_ab << endl;

    H.getIJ()({0,1},{0,1})({0,0}).writeRemoteData(ij_buf);
    H.getAI()({1,0},{0,1})({0,0}).writeRemoteData(ai_buf);
    H.getIA()({0,1},{1,0})({0,0}).writeRemoteData(ia_buf);
    H.getAB()({1,0},{1,0})({0,0}).writeRemoteData(ab_buf);

    H.getIJ()({0,0},{0,0}) = H.getIJ()({0,1},{0,1});
    H.getAI()({0,0},{0,0}) = H.getAI()({1,0},{0,1});
    H.getIA()({0,0},{0,0}) = H.getIA()({0,1},{1,0});
    H.getAB()({0,0},{0,0}) = H.getAB()({1,0},{1,0});

    writeIntegrals(abcd_buf, VABCD);
    writeIntegrals(abci_buf, VABCI);
    writeIntegrals(abij_buf, VABIJ);
    writeIntegrals(aibj_buf, VAIBJ);
    writeIntegrals(aijk_buf, VAIJK);
    writeIntegrals(ijkl_buf, VIJKL);

    H.getABCD()({2,0},{2,0})["ABCD"]  = 0.5*H.getABCD()({1,0},{1,0})["ABCD"];
    H.getABCD()({0,0},{0,0})["abcd"]  = 0.5*H.getABCD()({1,0},{1,0})["abcd"];

    H.getABCI()({2,0},{1,1})["ABCI"]  =     H.getABCI()({1,0},{1,0})["ABCI"];
    H.getABCI()({1,0},{0,1})["AbcI"]  =    -H.getABCI()({1,0},{1,0})["bAcI"];
    H.getABCI()({0,0},{0,0})["abci"]  =     H.getABCI()({1,0},{1,0})["abci"];

    H.getABIJ()({2,0},{0,2})["ABIJ"]  = 0.5*H.getABIJ()({1,0},{0,1})["ABIJ"];
    H.getABIJ()({0,0},{0,0})["abij"]  = 0.5*H.getABIJ()({1,0},{0,1})["abij"];

    H.getAIBJ()({1,1},{1,1})["AIBJ"]  =     H.getAIBJ()({1,0},{1,0})["AIBJ"];
    H.getAIBJ()({1,1},{1,1})["AIBJ"] -=     H.getABIJ()({1,0},{0,1})["ABJI"];
    H.getAIBJ()({0,1},{0,1})["aIbJ"]  =     H.getAIBJ()({1,0},{1,0})["aIbJ"];
    H.getAIBJ()({1,0},{0,1})["AibJ"]  =    -H.getABIJ()({1,0},{0,1})["AbJi"];
    H.getAIBJ()({0,1},{1,0})["aIBj"]  =    -H.getABIJ()({1,0},{0,1})["aBjI"];
    H.getAIBJ()({0,0},{0,0})["aibj"]  =     H.getAIBJ()({1,0},{1,0})["aibj"];
    H.getAIBJ()({0,0},{0,0})["aibj"] -=     H.getABIJ()({1,0},{0,1})["abji"];

    H.getAIJK()({1,1},{0,2})["AIJK"]  =     H.getAIJK()({1,0},{0,1})["AIJK"];
    H.getAIJK()({0,1},{0,1})["aIJk"]  =    -H.getAIJK()({1,0},{0,1})["aIkJ"];
    H.getAIJK()({0,0},{0,0})["aijk"]  =     H.getAIJK()({1,0},{0,1})["aijk"];

    H.getIJKL()({0,2},{0,2})["IJKL"]  = 0.5*H.getIJKL()({0,1},{0,1})["IJKL"];
    H.getIJKL()({0,0},{0,0})["ijkl"]  = 0.5*H.getIJKL()({0,1},{0,1})["ijkl"];

    H.getIJAK()({0,2},{1,1})["JKAI"]  =     H.getAIJK()({1,1},{0,2})["AIJK"];
    H.getIJAK()({0,1},{1,0})["JkAi"]  =     H.getAIJK()({1,0},{0,1})["AiJk"];
    H.getIJAK()({0,1},{0,1})["JkaI"]  =     H.getAIJK()({0,1},{0,1})["aIJk"];
    H.getIJAK()({0,0},{0,0})["jkai"]  =     H.getAIJK()({0,0},{0,0})["aijk"];

    H.getAIBC()({1,1},{2,0})["AIBC"]  =     H.getABCI()({2,0},{1,1})["BCAI"];
    H.getAIBC()({1,0},{1,0})["AiBc"]  =     H.getABCI()({1,0},{1,0})["BcAi"];
    H.getAIBC()({0,1},{1,0})["aIBc"]  =     H.getABCI()({1,0},{0,1})["BcaI"];
    H.getAIBC()({0,0},{0,0})["aibc"]  =     H.getABCI()({0,0},{0,0})["bcai"];

    H.getIJAB()({0,2},{2,0})["IJAB"]  =     H.getABIJ()({2,0},{0,2})["ABIJ"];
    H.getIJAB()({0,1},{1,0})["IjAb"]  =     H.getABIJ()({1,0},{0,1})["AbIj"];
    H.getIJAB()({0,0},{0,0})["ijab"]  =     H.getABIJ()({0,0},{0,0})["abij"];

    return true;
}

}
}

static const char* spec = R"!(

filename?
    string FCIDUMP

)!";

INSTANTIATE_SPECIALIZATIONS(aquarius::op::FCIDUMP);
REGISTER_TASK(aquarius::op::FCIDUMP<double>,"fcidump",spec);
