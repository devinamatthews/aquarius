#include "autocc.hpp"

namespace aquarius
{
namespace autocc
{

void compare(Diagram& d1, Diagram& d2)
{
    int width = 70;
    bool ok = true;

    //cout << "1:\n" << d1 << "2:\n" << d2;
    //cout << d1.numTerms() << " " << d2.numTerms() << endl;

    for (int i = 0;i < min(d1.numTerms(), d2.numTerms());i++)
    {
        ok = ok && d1[i] == d2[i] && d1[i].getFactor() == d2[i].getFactor();
        string s1, s2;
        ostringstream oss;
        oss << d1[i];
        s1 = oss.str();
        oss.str("");
        oss << d2[i];
        s2 = oss.str();
        cout << string(width-s1.size(), ' ') << s1 << string(width-s2.size(), ' ') << s2 << endl;
    }

    if (d1.numTerms() > d2.numTerms())
    {
        ok = false;
        for (int i = d2.numTerms();i < d1.numTerms();i++)
        {
            string s;
            ostringstream oss;
            oss << d1[i];
            s = oss.str();
            cout << string(width-s.size(), ' ') << s << endl;
        }
    }
    else if (d1.numTerms() < d2.numTerms())
    {
        ok = false;
        for (int i = d1.numTerms();i < d2.numTerms();i++)
        {
            string s;
            ostringstream oss;
            oss << d2[i];
            s = oss.str();
            cout << string(2*width-s.size(), ' ') << s << endl;
        }
    }

    if (!ok)
    {
        Diagram d3 = d1 - d2;
        cout << endl << "Diagrams are not equal" << endl << endl;
        for (int i = 0;i < d3.numTerms();i++)
        {
            string s;
            ostringstream oss;
            oss << d3[i];
            s = oss.str();
            cout << string(3*width/2-s.size(), ' ') << s << endl;
        }
    }
    else
    {
        cout << endl << "Diagrams are equal" << endl;
    }
}

int main(int argc, char **argv)
{
    ifstream fd(argv[1]);

    int start = 1;
    if (argc > 2)
    {
        istringstream(string(argv[2])) >> start;
    }

    vector<string> lines;
    string line;
    while (getline(fd, line)) lines.push_back(line);

    int d = 1;
    int lineno = start;
    while (lineno <= lines.size())
    {
        cout << "Checking diagram " << d << " (line " << lineno << ")" << endl << endl;

        //skip first line
        vector<string> diagram;
        while (true)
        {
            lineno++;
            string line = lines[lineno-1];
            if (line.find_first_not_of(" \t\n") == string::npos) break;
            diagram.push_back(line);
            cout << line << endl;
        }
        cout << endl;
        lineno++;

        cout << "Spin-orbital:" << endl;
        string so = lines[lineno-1];
        lineno++;
        int n = 0;
        if (so.find_first_of("AaIi") != string::npos) n++;
        if (so.find_first_of("BbJj") != string::npos) n++;
        if (so.find_first_of("CcKk") != string::npos) n++;
        if (so.find_first_of("DdLl") != string::npos) n++;
        cout << so << endl << endl;

        cout << "RHF:" << endl;
        vector<string> rhf;
        while (true)
        {
            if (lines[lineno-1].find('_') == string::npos) lineno++;
            string line = lines[lineno-1];
            if (line.find_first_not_of(" \t\n") == string::npos) break;
            if (line.find('_') != string::npos)
            {
                string mask(line);
                lineno++;
                line = lines[lineno-1];
                size_t i = string::npos;
                while ((i = mask.find('_', i+1)) != string::npos)
                {
                    line[i] = toupper(line[i]);
                }
            }
            rhf.push_back(line);
            cout << line << endl;
        }
        cout << endl;

        cout << "UHF:" << endl;
        vector<vector<string>> uhf(n+1);
        for (int k = 0;k <= n;k++)
        {
            while (true)
            {
                lineno++;
                string line = lines[lineno-1];
                if (line.find_first_not_of(" \t\n") == string::npos) break;
                uhf[k].push_back(line);
                cout << line << endl;
                if (lineno >= lines.size()) break;
            }
            cout << endl;
        }

        cout << "Checking RHF: ";
        string lower = string("aibjckdl").substr(0, 2*((n+1)/2));
        string upper = string("AIBJCKDL").substr(0, 2*((n+1)/2));
        string newso(so);
        Diagram d1(Diagram::SPINORBITAL, vector<string>(1, translate(newso, lower, upper)));
        d1.convert(Diagram::SKELETON).fixorder();
        for (vector<string>::iterator it = rhf.begin();it != rhf.end();++it) translate(*it, lower, upper);
        Diagram d2(Diagram::RHF, rhf);
        d2.convert(Diagram::SKELETON).fixorder();
        //cout << "1:\n" << d1 << endl;
        //cout << "2:\n" << d2 << endl;
        Diagram diff = d1 - d2;
        if (diff.numTerms() > 0)
        {
            cout << "failed" << endl;
            compare(d1, d2);
            return 1;
        }
        else
        {
            cout << "passed" << endl;
        }
        cout << endl;

        for (int k = 0;k <= n;k++)
        {
            cout << "Checking UHF " << string(n-k, 'A') << string(k, 'B') << ": ";
            string lower = string("aibjckdl").substr(0, 2*(n-k));
            string upper = string("AIBJCKDL").substr(0, 2*(n-k));
            //cout << "\nlower: " << lower << endl;
            //cout << "upper: " << upper << endl;
            string newso(so);
            Diagram d1(Diagram::SPINORBITAL, vector<string>(1, translate(newso, lower, upper)));
            d1.convert(Diagram::UHF).fixorder();
            Diagram d2(Diagram::UHF, uhf[k]);
            d2.fixorder();
            Diagram diff = d1 - d2;
            if (diff.numTerms() > 0)
            {
                cout << "failed" << endl;
                compare(d1, d2);
                return 1;
            }
            else
            {
                cout << "passed" << endl;
            }
        }
        cout << endl;

        //skip last line
        //skip -------...--------
        lineno += 2;
        if (lineno >= lines.size()) break;
        d++;
    }

    fd.close();

    return 0;
}

}
}
