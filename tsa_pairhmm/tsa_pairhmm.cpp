/*
Copyright (C) 2020 EMBL - European Bioinformatics Institute
Contact: goldman@ebi.ac.uk, cwalker@ebi.ac.uk

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*
This file is modified from the original, non-probabilistic implementation of
this method, which is distributed under the GNU General Public License v3.0
and is available at: https://github.com/ariloytynoja/fpa
*/

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

class TSApHMM
{

    /************************* global variables ******************************/

    vector<int> index1;
    vector<int> index2;
    vector<int> rindex1;
    vector<int> rindex2;

    vector<int> fseq1;
    vector<int> fseq2;

    vector<int> seq1;
    vector<int> rev1;
    vector<int> seq2;

    vector<int> mask1;

    int slg;
    int fsl1;
    int fsl2;


    int sl1;
    int sl2;
    int start1;
    int end1;
    int start2;
    int end2;
    int true_start2;
    int true_end2;

    int clus_start1;
    int clus_end1;
    int clus_start2;
    int clus_end2;

    string chrom;
    int chrom_start;

    string qry_name;
    string ref_name;

    //mod Store the scores for the forward/ts alignments
    float fwd_sco;
    float ts_sco;
    float traceback_max_sco;

    //mod Store the forward alignment sequence lengths
    int fwd_qry_seq_len;
    int fwd_ref_seq_len;

    //mod Store the template switch alignment sequence lengths
    int ts_qry_seq_len;
    int ts_ref_seq_len;

    struct Fasta_entry
    {
        string name;
        string sequence;
        int length;
    };

    struct switchPoint
    {
        int i;
        int j;
    };

    struct seqCoordinate
    {
        int pos_x;
        int pos_y;
        int matrix;
    };

    template<typename T>
    struct Array2D
    {
        private:
            int width;
            int org_width;
            int org_height;
            T * data;
        public:
            T& operator() (int x, int y) { return data[y*width + x]; }
            Array2D(const int w, const int h) : width(w), org_width(w), org_height(h) { data = new T[w*h]; }
            void resize(int nw, int nh) { if(nw*nh <= org_width*org_height) { width=nw; } else { delete [] data; data = new T[nw*nh]; org_width = nw; org_height = nh; } }
            ~Array2D() { delete [] data; }
    };

    enum Move_ptr {match=-1, xgap=-2, ygap=-3, none=-4};

    double maximum( float a, float b, float c ) { return max( max(a,b), c ) ; }
    
    /************************* global variables ******************************/



    /******************** command-line argument ******************************/

    po::options_description full_desc;
    po::options_description desc;
    po::options_description desc_dev;
    po::variables_map vm;

    bool pair_data;
    bool maximise_score;
    bool maximise_length;
    bool maximise_length_alt;
    bool verbose;

    void read_command_line_arguments(int argc, char *argv[])
    {

        desc.add_options()
                ("ref",  po::value<string>(), "reference (FASTA)")
                ("qry",  po::value<string>(), "query (FASTA)")
                ("pair", po::value<string>(), "sequence pair (FASTA)")
                ("qry-start", po::value<int>(), "query start")
                ("qry-end",   po::value<int>(), "query end")
                ("ref-start", po::value<int>(), "reference start")
                ("ref-end",   po::value<int>(), "reference end")
                ("ref-flank",   po::value<int>()->default_value(100), "reference flanking region")
                ("pair-start", po::value<int>(), "alignment start")
                ("pair-end",   po::value<int>(), "alignment end")
                ("scan","scan for mutation hotspots")
                ("scan-start", po::value<int>(), "scan start")
                ("scan-end",   po::value<int>(), "scan end")
                ("scan-flank",   po::value<int>()->default_value(40), "scan flanking region")
                ("scan-window-width",   po::value<int>()->default_value(10), "scan window width")
                ("scan-window-limit",   po::value<int>()->default_value(2), "minimum differences in scan window")
                ("switch-flank",   po::value<int>()->default_value(40), "switch event flanks")
                ("reverse","reverse sequence (control)")
                ("perfect-copy","copy region has to be identical")
                ("max-event-length",   po::value<int>()->default_value(1000), "maximum event length")
                ("cluster-annotation",  po::value<string>(), "output alignment with cluster annotation")
                ("print-file", po::value<string>(), "coordinates of events to print")
                ("long-output","long output")
                ("ref-name",  po::value<string>(), "reference name (for output)")
                ("qry-name",  po::value<string>(), "query name (for output)")
                ("aligner",  po::value<string>()->default_value("Input"), "aligner")
                ("swap-pair","swap sequence pair<")
                ("verbose", "verbose output")
                ("last-match", "the value to begin traceback from in the template switch alignment must be a match")
                ("ref_chrom", po::value<string>(), "name of species to get chromosome number from")
                ("rho", po::value<float>()->default_value(0.14), "rho value (indels relative to mismatches")
                ("divergence", po::value<float>()->default_value(0.01), "divergence time between the two species")
                ("lambda", po::value<float>()->default_value(20), "mean indel length, used to calculate epsilon")
                ("ts_n_switches", po::value<int>()->default_value(2750), "expected number of template switches genome-wide")
                ("ts_n_clusters", po::value<long>()->default_value(7937450), "number of mutation clusters genome-wide")
                ("ts_23_fragment_length", po::value<int>()->default_value(10), "estimate of the average 2->3 fragment length")

        ;
        desc_dev.add_options()
                ("max-23-score", "maximise score for 2-3 -fragment")
                // use to test alternate values of sigma
                ("max-23-length", "maximise length for 2-3 -fragment")
                ("max-23-length-alt", "maximise length for 2-3 -fragment")
                ("debug","debug")
                ("debug_scores","print score matrices as required")
                ("force-overlap","force overlap with cluster")
        ;

        full_desc.add(desc).add(desc_dev);
        po::store(po::command_line_parser(argc, argv).options(full_desc).run(), vm);

        if( ( not this->arg_is("ref") || not this->arg_is("qry") ) && not this->arg_is("pair") )
        {
            this->arg_error();
        }

        if(this->arg_is("verbose"))
            verbose= true;

        if(this->arg_is("max-23-score"))
            maximise_score = true;

        if(this->arg_is("max-23-length"))
            maximise_length = true;
        
        if(this->arg_is("max-23-length-alt"))
            maximise_length_alt = true;
    }


    void arg_error()
    {
        stringstream ss;
        ss<< "\nMany-to-One Bidirectional Alignment pairHMM\ncommand-line options\n" << desc << "\n";
        cout<<ss.str();
        exit(1);
    }

    bool arg_is(string name) { return vm.count(name); }

    const po::variable_value & arg_get(const string & name) { return vm[name]; }

    /******************** command-line argument ******************************/



    /********************  file input / output  ******************************/

    void read_data(string *datafile,Fasta_entry *entry)
    {
        ifstream input(datafile->c_str(), ios::in);
        if (!input) { return; }

        string temp, sequence = "";

        while(!input.eof())
        {
            getline(input, temp, '\n');

            if(temp[0] == '>')
            {
                entry->name = temp;
                entry->name.erase(entry->name.begin());
            }
            else sequence += temp;
        }

        entry->sequence = sequence;
        entry->length = sequence.length();
        input.close();
    }

    void read_pair_data(string *datafile,Fasta_entry *first,Fasta_entry *second)
    {
        ifstream input(datafile->c_str(), ios::in);
        if (!input) { return; }

        string temp, sequence = "";
        bool first_seq = true;

        while(!input.eof())
        {
            getline(input, temp, '\n');

            if(temp[0] == '>')
            {
                if(first_seq)
                {
                    first->name = temp;
                    first->name.erase(first->name.begin());
                    first_seq = false;
                }
                else
                {
                    first->sequence = sequence;
                    first->length = sequence.length();
                    sequence = "";
                    second->name = temp;
                    second->name.erase(second->name.begin());
                }
            }
            else sequence += temp;
        }

        second->sequence = sequence;
        second->length = sequence.length();

        if(arg_is("swap-pair"))
        {
            Fasta_entry tmp;
            tmp.name = second->name;
            tmp.sequence = second->sequence;

            second->name = first->name;
            second->sequence = first->sequence;
            second->length = second->sequence.length();

            first->name = tmp.name;
            first->sequence = tmp.sequence;
            first->length = first->sequence.length();

        }
        input.close();
    }



    /********************  file input / output  ******************************/



    /********************    alignment stuff    ******************************/

    void build_indeces(string *s1,string *s2)
    {

        if(s1->length() != s2->length())
        {
            cout<<"expecting aligned sequences. exiting.\n\n";
            exit(0);
        }

        slg = s1->length();

        index1.reserve(slg);
        index2.reserve(slg);
        rindex1.reserve(slg);
        rindex2.reserve(slg);

        fseq1.reserve(slg);
        fseq2.reserve(slg);

        int p1=0; int p2=0;

        for(int i=0;i<slg;i++)
        {
            index1.push_back(p1);
            index2.push_back(p2);

            if(s1->at(i) != '-')
            {
                rindex1.push_back(i);
                p1++;
            }
            if(s2->at(i) != '-')
            {
                rindex2.push_back(i);
                p2++;
            }
        }

        fsl1 = p1;
        fsl2 = p2;

        rindex1.resize(p1);
        rindex2.resize(p2);


        string alpha = "ACGTN";

        int ci;

        string::iterator sit = s1->begin();
        for(;sit != s1->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<5)
                fseq1.push_back(ci);
            else
                fseq1.push_back(-1);
        }

        sit = s2->begin();
        for(;sit != s2->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<5)
                fseq2.push_back(ci);
            else
                fseq2.push_back(-1);
        }

    }

    void build_sequences(string *s1,string *s2)
    {

        seq1.reserve(s1->length());
        rev1.reserve(s1->length());
        seq2.reserve(s2->length());
        mask1.reserve(s1->length());

        string alpha = "ACGTN";

        int ci;
        int p1=0; int p2=0;

        string::iterator sit = s1->begin();
        string::iterator sit2 = s2->begin();
        for(;sit != s1->end() && sit2 != s2->end();sit++,sit2++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<5)
            {
                seq1.push_back(ci);

                int m = 0;
                if(islower(*sit))
                    m += 1;
                if(islower(*sit2))
                    m += 2;

                mask1.push_back(m);

                p1++;
            }
        }


        sit = s1->begin();
        if(arg_is("reverse"))
        {
            for(;sit != s1->end();sit++)
            {
                ci = alpha.find(toupper(*sit));
                if (ci>=0 && ci<5)
                    rev1.push_back(ci);

            }
        }
        else
        {
            for(;sit != s1->end();sit++)
            {
                ci = alpha.find(toupper(*sit));
                if (ci>=0 && ci<4)
                    rev1.push_back(3-ci);
                else if (ci==4)
                    rev1.push_back(ci);

            }
        }

        sit = s2->begin();
        for(;sit != s2->end();sit++)
        {
            ci = alpha.find(toupper(*sit));
            if (ci>=0 && ci<5)
            {
                seq2.push_back(ci);
                p2++;
            }
        }

        seq1.resize(p1);
        rev1.resize(p1);
        seq2.resize(p2);
        mask1.resize(p1);


    }


    void set_two_fragments()
    {
        start1 = 0;
        end1 = seq1.size();

        if(arg_is("qry-start"))
        {
            int s = arg_get("qry-start").as<int>();
            if(s>=0 && s<end1)
                start1 = s;
        }
        if(arg_is("qry-end"))
        {
            int s = arg_get("qry-end").as<int>();
            if(s>start1 && s<end1)
                end1 = s;
        }

        sl1 = end1-start1;


        start2 = 0;
        end2 = seq2.size();

        if(arg_is("ref-start"))
        {
            int s = arg_get("ref-start").as<int>();
            if(s>=0 && s<end2)
                start2 = s;
        }
        if(arg_is("ref-end"))
        {
            int s = arg_get("ref-end").as<int>();
            if(s>start2 && s<end2)
                end2 = s;
        }

        int flank = arg_get("ref-flank").as<int>();

        true_start2 = start2;
        true_end2 = end2;

        if(flank>0 && start2-flank>=0)
            start2 = start2-flank;
        if(flank>0 && end2+flank<(int)seq2.size())
            end2 = end2+flank;

        sl2 = end2-start2;

    }

    void set_alignment_fragments()
    {

        int ps = 0;
        int pe = slg-1;

        if(arg_is("pair-start"))
        {
            int s = arg_get("pair-start").as<int>();
            if(s>=0 && s<slg)
                ps = s;
        }
        if(arg_is("pair-end"))
        {
            int s = arg_get("pair-end").as<int>();
            if(s>ps && s<slg)
                pe = s;
        }

        start1 = index1.at(ps);
        end1 = index1.at(pe);

        sl1 = end1-start1;

        int flank = 100;

        if(arg_is("ref-flank"))
        {
            flank = arg_get("ref-flank").as<int>();
        }

        start2 = index2.at(ps)-flank;
        end2 = index2.at(pe)+flank;

        if(start2 < 0)
            start2 = 0;

        if(end2>(int)seq2.size())
            end2 = seq2.size();

        sl2 = end2-start2;
    }


    float substitution_score_prob(float i, float j)
    {
        //cw match,mismatch emission probabilities
        float divergence = arg_get("divergence").as<float>();
        float log_emit_match    = log((0.75 * exp(-divergence)) + 0.25);
        float log_emit_mismatch = log(0.25 - 0.25 * exp(-divergence));
        
        if(seq1.at(i-1)<4 && seq1.at(i-1) == seq2.at(j-1))
            return log_emit_match;
        else
            return log_emit_mismatch;
    }

    
    float rev_substitution_score_prob(float i, float j)
    {
        float divergence = arg_get("divergence").as<float>();
        float log_emit_match    = log((0.75 * exp(-divergence)) + 0.25);
        float log_emit_mismatch = log(0.25 - 0.25 * exp(-divergence));

        if(arg_is("perfect-copy"))
        {
            if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                return log_emit_match;
            else
                return -100000;
        }
        else
        {
            if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                return log_emit_match;
            else
                return log_emit_mismatch;
        }
    }

    float substitution_score(float i, float j)
    {
        if(seq1.at(i-1)<4 && seq1.at(i-1) == seq2.at(j-1))
            return 1;
        else
            return -1;
    }


    int rev_substitution_score(int i, int j)
    {
        if(arg_is("perfect-copy"))
        {
            if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                return 1;
            else
                return -10000;
        }
        else
        {
            if(rev1.at(i-1)<4 && rev1.at(i-1) == seq2.at(j-1))
                return 1;
            else
                return -1;
        }
    }



    /********************    alignment stuff    ******************************/


    /********************    alignment output   ******************************/

    void print_switch_process(vector<seqCoordinate> *path,vector<switchPoint> *points)
    {
        string alpha = "ACGTN";

        string out1("");
        string out2("");
        string out2_gaps;
        string out3(" ");


        int point2 = points->at(1).j;
        int point3 = points->at(2).j;

        int p=path->size()-1;
        int site1 = path->at(p).pos_x;
        int site2 = path->at(p).pos_y;
        int site_mat = path->at(p).matrix;

        int first_site2 = site2;
        if(first_site2<0)
        {
            for(int i=p-1;i>0 && first_site2<0;i--)
                first_site2 = path->at(i).pos_y;
        }


        string qry(" ");
        string ref(" ");
        string rref(" ");

        int ref_pos = site2+1;
        int ref_end = 0;
        if(point3<site2)
        {  
            qry += " ";
            ref += " ";
            rref += " ";
            out1 += " ";
            out3 += " ";

            for(int i=point3;i<site2;i++)
            {
                out1 += " ";
                out3 += " ";
                ref += alpha.at(seq2.at(i-1));
                if(seq2.at(i-1)<4)
                    rref += alpha.at(3-seq2.at(i-1));
                else
                    rref += "N";
                ref_end = i-1;
            }
        }

        out1 += "\b L ";

        while(site_mat==1)
        {
            if(site1>=0 && site1<=(int)seq1.size())
            {
                out1+=alpha.at(seq1.at(site1-1));
                qry+=alpha.at(seq1.at(site1-1));
            }
            else
            {
                out1+="-";
                qry+="-";
            }
            if(site2>=0 && site2<=(int)seq2.size())
            {
                ref+=alpha.at(seq2.at(site2-1));
                if(seq2.at(site2-1)<4)
                    rref+=alpha.at(3-seq2.at(site2-1));
                else
                    rref += "N";
                ref_pos = site2+1;
            }
            else
            {
                ref+="-";
                rref+="-";
                out2_gaps+=" ";
                out3+=" ";
            }
            ref_end = site2-1;

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
            site_mat = path->at(p).matrix;
        }
        out1+=string(" 1");

        int first_site3 = site2;
        if(first_site3<0)
        {
            for(int i=p-1;i>0 && first_site3<0;i--)
                first_site3 = path->at(i).pos_y;
        }

        qry+="1 3";

        int last_y=0;

        while(site_mat==2)
        {
            if(site1>=0 && site1<=(int)seq1.size())
            {
                out2=alpha.at(seq1.at(site1-1))+out2;
                qry+=alpha.at(seq1.at(site1-1));
            }
            else
            {
                cout<<"error!\n";
            }

            last_y = site2;

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
            site_mat = path->at(p).matrix;
        }


        if(last_y>1)
            out2_gaps+=" ";

        if(ref_pos<=0)
            ref_pos=1;
        for(int i=ref_pos;i<site2 && i<(int)seq2.size();i++)
        {
            ref+=alpha.at(seq2.at(i-1));
            if(seq2.at(i-1)<4)
                rref += alpha.at(3-seq2.at(i-1));
            else
                rref += "N";
            ref_end = i-1;
        }

        out2=out2_gaps+string("\b3 ")+out2+string(" 2");
        for(int i=first_site2;i<last_y-1;i++)
        {
            out2=string(" ")+out2;
        }


        for(int i=first_site2;i<site2-1;i++)
        {
            out3+=string(" ");
            if(seq2.at(i)<0)
                out3+=string(" ");
        }

        out3+=string("\b4 ");

        qry+="2 4";

        while(site_mat==3 )
        {
            if(site1>=0 && site1<=(int)seq1.size())
            {
                out3+=alpha.at(seq1.at(site1-1));
                qry+=alpha.at(seq1.at(site1-1));
            }
            else
            {
                out3+="-";
                qry+="-";
            }
            if(site2>=0 && site2<=(int)seq2.size() && site2-1>ref_end)
            {
                ref+=alpha.at(seq2.at(site2-1));
                if(seq2.at(site2-1)<4)
                    rref+=alpha.at(3-seq2.at(site2-1));
                else
                    rref += "N";
            }
            else if( site2-1>ref_end )
            {
                ref+="-";
                rref+="-";
            }

            p--;
            if(p>=0)
            {
                site1 = path->at(p).pos_x;
                site2 = path->at(p).pos_y;
                site_mat = path->at(p).matrix;
            }
            else
                break;
        }
        site2++;
        if(point2>site2)
        {
            for(;site2<=point2;site2++)
            {
                ref+=alpha.at(seq2.at(site2-1));
                if(seq2.at(site2-1)<4)
                    rref+=alpha.at(3-seq2.at(site2-1));
                else
                    rref += "N";
            }
        }
        out3+=string(" R");

        string qry0;
        for(int i=start1;i<end1;i++)
            qry0+=alpha.at(seq1.at(i));

        string ref0;
        for(int i=start2;i<end2;i++)
            ref0+=alpha.at(seq2.at(i));

        string ref00;
        for(int i=true_start2;i<true_end2;i++)
            ref00+=alpha.at(seq2.at(i));

        if(!arg_is("scan"))
        {
            cout<<endl<<"chr"<<chrom<<":"<<chrom_start+points->at(0).i+1<<"-"<<chrom_start+points->at(0).i-points->at(2).j+points->at(1).j+1<<endl<<endl;
        }

        if(verbose)
        {
            cout<<"Template switch process:"<<"\nF1: "<<out1<<"\nF3:  "<<out3<<"\nRF:  "<<ref<<"\nRR:  "<<rref<<"\nF2:  "<<out2<<"\n\n";
        }
    }



    void print_inversion_fragment(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path)
    {

        string alpha = "-ACGTN";
        string lowalpha = "-acgtn";
        int flanking = arg_get("switch-flank").as<int>();

        bool reverse = arg_is("reverse");

        int start_i = max(points->at(0).i-flanking,0);
        int stop_i = points->at(0).i;

        int p=path->size()-1;

        int site1 = path->at(p).pos_x;
        int site2 = path->at(p).pos_y;

        while(site1<=start_i)
        {
            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }

        vector<int> seq1_frag1;
        vector<int> seq2_frag1;
        vector<bool> low_frag1;

        while(site1<=stop_i)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
               seq1.at(site1-1) != seq2.at(site2-1) )
                low_frag1.push_back(true);
            else if(site2<0)
                low_frag1.push_back(true);
            else
                low_frag1.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
                seq1_frag1.push_back(seq1.at(site1-1));
            else
                seq1_frag1.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_frag1.push_back(seq2.at(site2-1));
            else
                seq2_frag1.push_back(-1);

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }

        vector<int> seq1_frag2;
        vector<int> seq2_frag2;
        vector<bool> low_frag2;

        stop_i = points->at(2).i;

        int sA=0; int sC=0; int sG=0; int sT=0;
        while(site1<=stop_i)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
              ( (reverse && seq1.at(site1-1) != seq2.at(site2-1) ) ||
                ( not reverse && seq1.at(site1-1) != 3-seq2.at(site2-1) ) ) )
                low_frag2.push_back(true);
            else if(site2<0)
                low_frag2.push_back(true);
            else
                low_frag2.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
            {
                seq1_frag2.push_back(seq1.at(site1-1));
                int c = seq1.at(site1-1);
                if(c==0) sA=1;
                if(c==1) sC=1;
                if(c==2) sG=1;
                if(c==3) sT=1;
            }
            else
                seq1_frag2.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_frag2.push_back(3-seq2.at(site2-1));
            else
                seq2_frag2.push_back(-1);

            p--;
            site1 = path->at(p).pos_x;
            site2 = path->at(p).pos_y;
        }
        // base count in 2->3 region
        int sumNuc = sA+sC+sG+sT;

        vector<int> seq1_frag3;
        vector<int> seq2_frag3;
        vector<bool> low_frag3;

        stop_i = min(points->at(3).i+flanking,(int)seq1.size());

        while(site1<stop_i && p>=0)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
               seq1.at(site1-1) != seq2.at(site2-1) )
                low_frag3.push_back(true);
            else if(site2<0)
                low_frag3.push_back(true);
            else
                low_frag3.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
                seq1_frag3.push_back(seq1.at(site1-1));
            else
                seq1_frag3.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_frag3.push_back(seq2.at(site2-1));
            else
                seq2_frag3.push_back(-1);

            p--;
            if(p>=0)
            {
                site1 = path->at(p).pos_x;
                site2 = path->at(p).pos_y;
            }
        }

        vector<int> seq1_fwd;
        vector<int> seq2_fwd;
        vector<bool> low_fwd;

        p=fwd_path->size()-1;

        site1 = fwd_path->at(p).pos_x;
        site2 = fwd_path->at(p).pos_y;

        while(site1<=start_i)
        {
            p--;
            site1 = fwd_path->at(p).pos_x;
            site2 = fwd_path->at(p).pos_y;
        }
        stop_i = min(points->at(3).i+flanking,(int)seq1.size());
        int epo_start1 = site1-1;
        int epo_stop1 = stop_i-1;

        while(site1<stop_i && p>=0)
        {
            if(site1>0 && site1<=(int)seq1.size() &&
               site2>0 && site2<=(int)seq2.size() &&
               seq1.at(site1-1) != seq2.at(site2-1) )
                low_fwd.push_back(true);
            else if(site2<0)
                low_fwd.push_back(true);
            else
                low_fwd.push_back(false);

            if(site1>0 && site1<=(int)seq1.size())
                seq1_fwd.push_back(seq1.at(site1-1));
            else
                seq1_fwd.push_back(-1);

            if(site2>0 && site2<=(int)seq2.size())
                seq2_fwd.push_back(seq2.at(site2-1));
            else
                seq2_fwd.push_back(-1);

            p--;
            if(p>=0)
            {
                site1 = fwd_path->at(p).pos_x;
                site2 = fwd_path->at(p).pos_y;
            }
        }


        float up_ident;
        float repeat_ident;
        float down_ident;
        float inv_ident;
        float fwd_ident;

        int inv_sum_length = 0;
        int inv_sum_ins = 0;
        int inv_sum_del = 0;
        int inv_sum_mis = 0;

        int sum1=0;
        for(int i=0;i<(int)seq1_frag1.size();i++)
        {
            if(seq1_frag1.at(i)==seq2_frag1.at(i))
                sum1++;
            else if(seq1_frag1.at(i)<0)
                inv_sum_del++;
            else if(seq2_frag1.at(i)<0)
                inv_sum_ins++;
            else
                inv_sum_mis++;

            inv_sum_length++;
         }

        if(seq1_frag1.size()>0)
            up_ident = float(sum1)/int(seq1_frag1.size());
        else
            up_ident = 0;

        int sum2=0;
        int pState=-1;
        bool hasCG=false;
        bool hasGC=false;
        for(int i=0;i<(int)seq1_frag2.size();i++)
        {
            if(not reverse && seq1_frag2.at(i)==seq2_frag2.at(i))
                sum2++;
            else if(reverse && seq1_frag2.at(i)>=0 && seq1_frag2.at(i)<4 && seq1_frag2.at(i)==3-seq2_frag2.at(i))
                sum2++;
            else if(seq1_frag2.at(i)<0)
                inv_sum_del++;
            else if(seq2_frag2.at(i)<0)
                inv_sum_ins++;
            else
                inv_sum_mis++;

            inv_sum_length++;

            if(pState==1 && seq1_frag2.at(i)==2)
                hasCG=true;
            else if(pState==2 && seq1_frag2.at(i)==1)
                hasGC=true;

            pState=seq1_frag2.at(i);
        }
        int CpG=0;
        if(hasCG)
            CpG+=1;
        if(hasGC)
            CpG+=2;

        repeat_ident = float(sum2)/int(seq1_frag2.size());


        int sum3=0;
        for(int i=0;i<(int)seq1_frag3.size();i++)
        {
            if(seq1_frag3.at(i)==seq2_frag3.at(i))
                sum3++;
            else if(seq1_frag3.at(i)<0)
                inv_sum_del++;
            else if(seq2_frag3.at(i)<0)
                inv_sum_ins++;
            else
                inv_sum_mis++;

            inv_sum_length++;
         }

        if(seq1_frag3.size()>0)
            down_ident = float(sum3)/int(seq1_frag3.size());
        else
            down_ident = 0;

        inv_ident = float(sum1+sum2+sum3)/int(seq1_frag1.size()+seq1_frag2.size()+seq1_frag3.size());

        int sum_mis = 0;
        int sum_ins = 0;
        int sum_del = 0;

        int sum4=0;
        for(int i=0;i<(int)seq1_fwd.size();i++)
        {
            if(seq1_fwd.at(i)==seq2_fwd.at(i))
                sum4++;
            else
            {
                if(seq1_fwd.at(i)<0)
                    sum_del++;
                else if(seq2_fwd.at(i)<0)
                    sum_ins++;
                else
                    sum_mis++;
            }
        }
        fwd_ident = float(sum4)/int(seq1_fwd.size());


        float epo_ident = 0;
        int mask_state = 0;

        int clus_ins = 0;
        int clus_del = 0;
        int clus_mis = 0;

        if(arg_is("scan") || arg_is("print-file"))
        {
            int epo_sum = 0; int epo_length = 0;
            this->fwd_compare_sequences(&epo_sum,&epo_length,epo_start1,epo_stop1);
            epo_ident = float(epo_sum)/epo_length;

            int m_start = clus_start1;
            if(m_start>0)
                m_start--;
            if(m_start>0)
                m_start--;

            if(clus_end1>=fsl1)
                clus_end1 = fsl1-1;

            int m_end = clus_end1;
            if(m_end+1<fsl1)
                m_end++;
            if(m_end+1<fsl1)
                m_end++;

            for(int i=m_start;i<m_end;i++)
            {
                if(mask1.at(i))
                {
                    if(mask1.at(i)>mask_state)
                        mask_state = mask1.at(i);
                }
            }

            int clus_start = rindex1.at(clus_start1);
            int clus_end = rindex1.at(clus_end1);

            if(clus_start1>0)
                clus_start = rindex1.at(clus_start1-1);
            if(clus_end1<(int)rindex1.size()-1)
                clus_end = rindex1.at(clus_end1+1);

            for(int i=clus_start;i<clus_end;i++)
            {
                if(fseq1.at(i)<0)
                    clus_del++;
                else if(fseq2.at(i)<0)
                    clus_ins++;
                else if(fseq1.at(i) != fseq2.at(i))
                    clus_mis++;
            }
        }
        
        // fragment sizes of L->1, 2->3, 4->R for output
        int f1_size, f2_size, f3_size;
        f1_size = seq1_frag1.size();
        f2_size = seq1_frag2.size();
        f3_size = seq1_frag3.size();
        //ts_qry_seq_len = f1_size + f2_size + f3_size;

        cout<<setprecision(3);

        if(arg_is("scan") || arg_is("print-file"))
        {
            if(arg_is("scan") || arg_is("long-output"))
            {
                int c_start = rindex1.size();
                if(clus_start1<(int)rindex1.size())
                    c_start = rindex1.at(clus_start1);

                if(verbose)
                    cout<<"chrom,clus_start_chrom,clus_start_align,clust_start1,clust_end1,sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref,iden_up,ident_rep,ident_down,ident_inv,ident_fwd,ident_epo,masked,sum_ins,sum_del,sum_mis,sum_nuc,CpG,clus_ins,clus_del,clus_mis,fwd_score,ts_score_local,ts_score_global,ts_ref_seq_len,ts_qry_seq_len,frag_L1_size,frag_23_size,frag_4R_size\n";
                cout<<chrom<<","<<chrom_start+clus_start1<<","<<c_start<<","<<clus_start1<<","<<clus_end1<<","<<points->at(0).i<<","<<points->at(0).j<<","<<points->at(1).j<<","<<points->at(2).j<<","<<points->at(3).j<<","
                    <<up_ident<<","<<repeat_ident<<","<<down_ident<<","<<inv_ident<<","<<fwd_ident<<","<<epo_ident<<","<<mask_state<<","<<sum_ins-inv_sum_ins<<","<<sum_del-inv_sum_del<<","<<sum_mis-inv_sum_mis<<","<<sumNuc<<","
                    <<CpG<<","<<clus_ins<<","<<clus_del<<","<<clus_mis<<","<<fwd_sco<<","<<ts_sco<<","<<traceback_max_sco<<","<<ts_ref_seq_len<<","<<ts_qry_seq_len<<","<<f1_size<<","<<f2_size<<","<<f3_size<<"\n";
            }
        }
        else
        {
            if(verbose)
                cout<<"sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref,iden_up,ident_rep,ident_down,ident_inv,ident_fwd,sum_ins,sum_del,sum_mis,sum_nuc,CpG\n";
            cout<<points->at(0).i<<","<<points->at(0).j<<","<<points->at(1).j<<","<<points->at(2).j<<","<<points->at(3).j<<","
                <<up_ident<<","<<repeat_ident<<","<<down_ident<<","<<inv_ident<<","<<fwd_ident<<","<<sum_ins<<","<<sum_del<<","<<sum_mis<<","<<sumNuc<<","<<CpG<<"\n";
        }


        if(verbose)
        {
            stringstream qry;
            stringstream ref;
            for(int i=0;i<(int)seq1_frag1.size();i++)
            {
                if(low_frag1.at(i))
                    qry<<lowalpha.at(seq1_frag1.at(i)+1);
                else
                    qry<<alpha.at(seq1_frag1.at(i)+1);
                ref<<alpha.at(seq2_frag1.at(i)+1);
            }
            ref<<"|";
            qry<<"|";
            for(int i=0;i<(int)seq1_frag2.size();i++)
            {
                if(low_frag2.at(i))
                    qry<<lowalpha.at(seq1_frag2.at(i)+1);
                else
                    qry<<alpha.at(seq1_frag2.at(i)+1);
                ref<<alpha.at(seq2_frag2.at(i)+1);
            }
            ref<<"|";
            qry<<"|";
            for(int i=0;i<(int)seq1_frag3.size();i++)
            {
                if(low_frag3.at(i))
                    qry<<lowalpha.at(seq1_frag3.at(i)+1);
                else
                    qry<<alpha.at(seq1_frag3.at(i)+1);
                ref<<alpha.at(seq2_frag3.at(i)+1);
            }
            ref<<"\n";
            qry<<"\n";

            stringstream fwd_qry;
            stringstream fwd_ref;
            for(int i=0;i<(int)seq1_fwd.size();i++)
            {
                if(low_fwd.at(i))
                    fwd_qry<<lowalpha.at(seq1_fwd.at(i)+1);
                else
                    fwd_qry<<alpha.at(seq1_fwd.at(i)+1);
                fwd_ref<<alpha.at(seq2_fwd.at(i)+1);
            }
            fwd_ref<<"\n";
            fwd_qry<<"\n";
            

            //if(arg_is("long-output"))
            cout<<"\nUnidirectional alignment (log-probability: "<<fwd_sco<<")"<<"\n"<<qry_name<<" "<<fwd_qry.str()<<ref_name<<" "<<fwd_ref.str()<<endl;

            cout<<"Template switch alignment (log-probability: "<<ts_sco<<")"<<"\n"<<qry_name<<" "<<qry.str()<<ref_name<<" "<<ref.str()<<endl;

            
            if(arg_is("long-output"))
            {
                cout<<"Input alignment:\n";

                int start1_epo = max( fwd_path->at(fwd_path->size()-1).pos_x-1, points->at(0).i-flanking ) ;
                int ss1 = fwd_path->at(0).pos_x+1;
                if(ss1<0)ss1=(int)seq1.size();
                int stop1_epo = min( points->at(3).i+flanking,min(ss1,(int)seq1.size()-1 ));

                cout<<qry_name<<" ";
                for(int i=rindex1.at(start1_epo);i<rindex1.at(stop1_epo)-1;i++)
                    if(fseq1.at(i)==fseq2.at(i))
                       cout<<string("-ACGTN").at(fseq1.at(i)+1);
                    else
                       cout<<string("-acgtn").at(fseq1.at(i)+1);
                cout<<endl<<ref_name<<" ";

                for(int i=rindex1.at(start1_epo);i<rindex1.at(stop1_epo)-1;i++)
                    cout<<string("-ACGTN").at(fseq2.at(i)+1);
                cout<<endl<<endl;
            }

        }

    }

    /********************    alignment output   ******************************/



    /********************    alignment itself   ******************************/

    void fwd_compare_sequences(int *identical,int *length,int epo_start1,int epo_stop1)
    {
        *identical = 0;
        *length = 0;

        for(int i=rindex1.at(epo_start1);i<rindex1.at(epo_stop1);i++)
        {
            if(fseq1.at(i) == fseq2.at(i))
                (*identical)++;
            (*length)++;
        }
    }

    void fwd_align_sequences(vector<seqCoordinate> *path,vector<seqCoordinate> *inv_path)
    {
        int p = inv_path->size()-1;
        int si = inv_path->at(p).pos_x;
        while(si<0 && p>=0)
            si = inv_path->at(--p).pos_x;

        p = inv_path->size()-1;
        int sj = inv_path->at(p).pos_y;
        while(sj<0 && p>=0)
            sj = inv_path->at(--p).pos_y;
        p=0;
        int ei = inv_path->at(0).pos_x;
        while(ei<0 && p<(int)inv_path->size())
            ei = inv_path->at(++p).pos_x;
        p=0;
        int ej = inv_path->at(p).pos_y;
        while(ej<0 && p<(int)inv_path->size())
            ej = inv_path->at(++p).pos_y;

        int tmp = min(si,ei);
        ei = max(si,ei);
        si = tmp;
        tmp = min(sj,ej);
        ej = max(sj,ej);
        sj = tmp;

        if(arg_is("scan") || arg_is("print-file"))
        {
            si = start1;
            ei = end1;
            sj = true_start2;
            ej = true_end2;
        }

        //cw Added two additional scoring matrices to allow affine gap penalties
        Array2D<double> fwd_mat1(ei-si+1,ej-sj+1);
        Array2D<double> fwd_mat2(ei-si+1,ej-sj+1);
        Array2D<double> fwd_mat3(ei-si+1,ej-sj+1);

        //cw Added two additional pointer matrices to allow affine gap penalties
        Array2D<int> fwd_ptr1(ei-si+1,ej-sj+1);
        Array2D<int> fwd_ptr2(ei-si+1,ej-sj+1);
        Array2D<int> fwd_ptr3(ei-si+1,ej-sj+1);

        fwd_qry_seq_len=ei-si;
        fwd_ref_seq_len=ej-sj;

        float large_neg = -10000;
        
        // variables to calculate state log probabilities
        double time_divergence = arg_get("divergence").as<float>();
        double rho = arg_get("rho").as<float>();
        double lambda = arg_get("lambda").as<float>();
        double delta = 1.0 - exp(-time_divergence*(rho/2));
        double eps = 1.0 - (1/lambda);

        // emission probabilities
        double log_emit_indel = log(0.25);

        // state log probabilities
        double log_match_indel = log(delta);
        double log_indel_match = log((1-eps) * (1-(2*delta)));
        double log_indel_indel    = log((1-eps) * delta);
        double log_match_extend   = log(1-(2*delta));
        double log_indel_extend   = log(eps + ((1-eps)*delta));
            
        bool debug_score_matrix = arg_is("debug_scores");


        for(int j=0;j<=ej-sj;j++)
        {
            // fwd_mat1
            for(int i=0;i<=ei-si;i++)
            {
                float max_score = large_neg;
                
                if(i==0 && j==0) {
                    fwd_mat1(i,j) = 0;
                    fwd_mat2(i,j) = log_emit_indel;
                    fwd_mat3(i,j) = log_emit_indel;
                    
                    fwd_ptr1(i,j) = none;
                    fwd_ptr2(i,j) = none;
                    fwd_ptr3(i,j) = none;

                } else {
                    if(i>0 && j>0)
                    {   
                        float sub_score = substitution_score_prob(si+i,sj+j);
                        float match_score = fwd_mat1(i-1,j-1) + sub_score + log_match_extend;
                        float xgap_score =  fwd_mat2(i-1,j-1) + sub_score + log_indel_match;
                        float ygap_score =  fwd_mat3(i-1,j-1) + sub_score + log_indel_match;
                        max_score = maximum(match_score, xgap_score, ygap_score);
                        if(match_score == max_score) {
                            fwd_ptr1(i,j) = match;
                        } else if(xgap_score == max_score) {
                            fwd_ptr1(i,j) = xgap;
                        } else if(ygap_score == max_score) {
                            fwd_ptr1(i,j) = ygap;
                        } else {
                            cout<<"i>0 & j>0 none match max, exiting"<<endl;
                            exit(0);
                        }
                        fwd_mat1(i,j) = max_score;
                    } else {
                        fwd_mat1(i,j) = large_neg;
                        fwd_ptr1(i,j) = none;
                    }
                    
                    if(i>0)
                    {
                        float match_open = fwd_mat1(i-1,j) + log_emit_indel + log_match_indel;
                        float y_open     = fwd_mat3(i-1,j) + log_emit_indel + log_indel_indel;
                        float x_extend   = fwd_mat2(i-1,j) + log_emit_indel + log_indel_extend;
                        max_score = maximum(match_open, x_extend, y_open);
                        if(match_open == max_score) {
                            fwd_ptr2(i,j) = match;
                        } else if(x_extend == max_score) {
                            fwd_ptr2(i,j) = xgap;
                        } else if(y_open == max_score) {
                            fwd_ptr2(i,j) = ygap;
                        } else {
                            cout<<"i>0 none match max, exiting"<<endl;
                            exit(0);
                        }
                        fwd_mat2(i,j) = max_score;
                    } else {
                        fwd_mat2(i,j) = large_neg;
                        fwd_ptr2(i,j) = none;
                    }

                    if(j>0)
                    {
                        float match_open = fwd_mat1(i,j-1) + log_match_indel; // + log_emit_indel;
                        float x_open     = fwd_mat2(i,j-1) + log_indel_indel; //+ log_emit_indel;
                        float y_extend   = fwd_mat3(i,j-1) + log_indel_extend; //+ log_emit_indel;
                        max_score = maximum(match_open, y_extend, x_open);
                        if(match_open == max_score) {
                            fwd_ptr3(i,j) = match;
                        } else if(y_extend == max_score) {
                            fwd_ptr3(i,j) = ygap;
                        } else if(x_open == max_score) {
                            fwd_ptr3(i,j) = xgap;
                        } else {
                            cout<<"j>0 none match max, exiting"<<endl;
                            exit(0);
                        }
                        fwd_mat3(i,j) = max_score;
                    } else {
                        fwd_mat3(i,j) = large_neg;
                        fwd_ptr3(i,j) = none;
                    }
                }
            }
            
        }
            
        if(debug_score_matrix) {
            cout<<"FWD_M"<<endl;
            for(int i=0;i<=ei-si;i++) {
                for(int j=0;j<=ej-sj;j++) {
                    cout<<fwd_mat1(i,j)<<",";
               }
            cout<<endl;
            }
        }

        int i=ei-si;
        int j=ej-sj;
        int current = 0;
        Array2D<int> *current_ptr = &fwd_ptr1;
        
        float max_sco = maximum(fwd_mat1(i,j), fwd_mat2(i,j), fwd_mat3(i,j));
        if(max_sco == fwd_mat1(i,j)) {
            current_ptr = &fwd_ptr1;
            current = 1;
        } else if(max_sco == fwd_mat2(i,j)) {
            current_ptr = &fwd_ptr2;
            current = 2;
        } else if(max_sco == fwd_mat3(i,j)) {
            current_ptr = &fwd_ptr3;
            current = 3;
        } else {
            cout<<"fwd alignment max not correct"<<endl;
            exit(0);
        }
        // Was used to print the forward matrix for visualisation
        /*  
        for(int i_p=0;i_p<=i;i_p++)
            {
                for(int j_p=0;j_p<=j;j_p++)
                {
                    cout<<fwd_mat3(i_p,j_p)<<",";
                }
                cout<<endl;
            }
        */
        fwd_sco = max_sco;
 
        for(;i>=0 || j>=0;) {
            seqCoordinate c;
            c.matrix = 1;
            c.pos_x = si+i;
            c.pos_y = sj+j;
            
            if(current==1) {
                if((*current_ptr)(i,j)==match) {
                    path->push_back(c);
                    i--;j--;
                    current_ptr = &fwd_ptr1;
                    current = 1;
                
                } else if((*current_ptr)(i,j)==xgap) {
                    path->push_back(c);
                    i--;j--;
                    current_ptr = &fwd_ptr2;
                    current = 2;
                
                } else if((*current_ptr)(i,j)==ygap) {
                    path->push_back(c);
                    i--;j--;
                    current_ptr = &fwd_ptr3;
                    current = 3;
                
                } else if((*current_ptr)(i,j)==none) {
                    current = 5;
                    break;
                
                } else {
                    cout<<"current ptr not correct"<<endl;
                    cout<<(*current_ptr)(i,j)<<endl;
                    exit(0);
                }


            } else if(current==2) {
                if((*current_ptr)(i,j)==match) {
                    c.pos_y = -1;
                    path->push_back(c);
                    i--;
                    current_ptr = &fwd_ptr1;
                    current = 1;
                
                } else if((*current_ptr)(i,j)==xgap) {
                    c.pos_y = -1;
                    path->push_back(c);
                    i--;
                    current_ptr = &fwd_ptr2;
                    current = 2;
                
                } else if((*current_ptr)(i,j)==ygap) {
                    c.pos_y = -1;
                    path->push_back(c);
                    i--;
                    current_ptr = &fwd_ptr3;
                    current = 3;
                
                } else if((*current_ptr)(i,j)==none) {
                    current = 5;
                    break;
                
                } else {
                    cout<<"current ptr not correct:"<<endl;
                    cout<<(*current_ptr)(i,j)<<endl;
                    exit(0);
                }

            } else if(current==3) {
                if((*current_ptr)(i,j)==match) {
                    c.pos_x = -1;
                    path->push_back(c);
                    j--;
                    current_ptr = &fwd_ptr1;
                    current = 1;

                } else if((*current_ptr)(i,j)==xgap) {
                    c.pos_x = -1;
                    path->push_back(c);
                    j--;
                    current_ptr = &fwd_ptr2;
                    current = 2;

                } else if((*current_ptr)(i,j)==ygap) {
                    c.pos_x = -1;
                    path->push_back(c);
                    j--;
                    current_ptr = &fwd_ptr3;
                    current = 3;

                } else if((*current_ptr)(i,j)==none) {
                    current = 5;
                    break;

                } else {
                    cout<<"current ptr not correct"<<endl;
                    cout<<(*current_ptr)(i,j)<<endl;
                    exit(0);
                }

            } else if(current==5) {
                break;

            } else {
                cout<<"current not equal 1,2,3"<<endl;
                exit(0);
            }

        }
        
        }

        void align_sequences(vector<seqCoordinate> *path,vector<switchPoint> *points,bool local=true)
        {
            int maxtermgap = 0;
            Array2D<double> mat1(sl1+1,sl2+1);
            Array2D<double> mat1_2(sl1+1,sl2+1);
            Array2D<double> mat1_3(sl1+1,sl2+1);

            Array2D<double> mat2(sl1+1,sl2+1);
            //Array2D<double> mat2_2(sl1+1,sl2+1);
            //Array2D<double> mat2_3(sl1+1,sl2+1);

            Array2D<double> mat3(sl1+1,sl2+1);
            Array2D<double> mat3_2(sl1+1,sl2+1);
            Array2D<double> mat3_3(sl1+1,sl2+1);

            Array2D<int> ts_ptr1(sl1+1,sl2+1);
            Array2D<int> ts_ptr1_2(sl1+1,sl2+1);
            Array2D<int> ts_ptr1_3(sl1+1,sl2+1);
            
            Array2D<int> ptr2(sl1+1,sl2+1);
            
            Array2D<int> ptr3(sl1+1,sl2+1);
            Array2D<int> ptr3_2(sl1+1,sl2+1);
            Array2D<int> ptr3_3(sl1+1,sl2+1);

            Array2D<double> sco2(sl1+1,sl2+1);
            Array2D<int> len2(sl1+1,sl2+1);
            Array2D<int> dist2(sl1+1,sl2+1);
            Array2D<int> jump_mat(sl1+1,sl2+1);

            Array2D<double> sco3(sl1+1,sl2+1);
            Array2D<int> len3(sl1+1,sl2+1);

            int large_neg = -100000;

            ts_qry_seq_len=sl1;
            ts_ref_seq_len=sl2;

            double time_divergence = arg_get("divergence").as<float>();
            double rho = arg_get("rho").as<float>();
            double lambda = arg_get("lambda").as<float>();
            double delta = 1.0 - exp(-time_divergence*(rho/2));
            double eps = 1.0 - (1/lambda);
            
            
            // template switch penalties
            double n_switches = arg_get("ts_n_switches").as<int>();
            double n_clusters = arg_get("ts_n_clusters").as<long>();
            double average_two_three_length = arg_get("ts_23_fragment_length").as<int>();

            double theta = n_switches / (n_clusters * ts_qry_seq_len);
            double sigma = 1 / average_two_three_length;
            
            // set m2 depending on preferring long 2->3 regions (for testing parameters)
            if(maximise_length) {
                sigma = delta;
            }
            
            if(maximise_length_alt) {
                sigma = 2 * delta;
            }

            // emission probabilities
            double log_emit_indel     = log(0.25);
            
            // state log probabilities
            double m1_to_i1d1  = log(delta);
            double i1d1_to_m1  = log((1-eps) * (1-(2*delta)-theta));
            double i1d1_to_m2  = log( (1-eps) * theta );
            double m1_extend   = log(1 - (2*delta) - theta);
            double m1_to_m2    = log(theta);

            double m2_extend   = log(1-sigma);

            double m2_to_m3    = log( sigma * (1-(2*delta)) );
            double m2_to_i3d3  = log( sigma * delta );

            double i3d3_to_m3 = log( (1-eps) * (1-(2*delta)) );
            double m3_extend = log( 1 - (2*delta) );
            double m3_to_i3d3 = log(delta);

            // indel parameters which are the same in i1,d1,i3,d3
            double indel_to_indel     = log((1-eps) * delta);
            double indel_extend   = log(eps + ((1-eps) * delta));

            // check if debug flags were specified
            bool debug_matrix = arg_is("debug");
            bool last_match = arg_is("last-match");
            bool debug_score_matrix = arg_is("debug_scores");


            if(debug_matrix)
            {   
                cout<<start1<<" "<<end1<<"; "<<start2<<" "<<end2<<"; "<<true_start2<<" "<<true_end2<<" | "<<clus_start1<<" "<<clus_end1<<"; "<<clus_start2<<" "<<clus_end2<<"\n";
                for(int i=0;i<=sl1;i++)
                     cout<<""<<string("-ACGTN").at(seq1.at(start1+i)+1);
                cout<<endl;
            }

            bool force_overlap = arg_is("force-overlap");

            float maxs = 0;
            int maxi = 0;
            int maxj = 0;

            for(int j=0;j<=sl2;j++)
            {
                if(debug_matrix) cout<<start2+j<<" "<<string("-ACGTN").at(seq2.at(start2+j)+1);
                // mat1
                for(int i=0;i<=sl1;i++)
                {
                    float max_score = large_neg;
                    float top_score = large_neg;

                    if(i==0 && j==0) {
                        ts_ptr1(0,0) = none;
                        ts_ptr1_2(0,0) = none;
                        ts_ptr1_3(0,0) = none;
                        
                        mat1(0,0) = 0;
                        mat1_2(0,0) = log_emit_indel + m1_to_i1d1;
                        mat1_3(0,0) = log_emit_indel + m1_to_i1d1;

                    } else if( (i==0 && j<=maxtermgap) || (j==0 && i<=maxtermgap) ) {
                        ts_ptr1(i,j) = none;
                        ts_ptr1_2(i,j) = none;
                        ts_ptr1_3(i,j) = none;
                        
                        mat1(i,j) = 0;
                        mat1_2(i,j) = log_emit_indel + m1_to_i1d1;
                        mat1_3(i,j) = log_emit_indel + m1_to_i1d1;
                    
                    } else if(not local && i==0 && j+start2==true_start2 ) { 
                        ts_ptr1(i,j) = none;
                        ts_ptr1_2(i,j) = none;
                        ts_ptr1_3(i,j) = none;
                        
                        mat1(i,j) = 0;
                        mat1_2(i,j) = log_emit_indel + m1_to_i1d1;
                        mat1_3(i,j) = log_emit_indel + m1_to_i1d1;

                    } else {
                        if(i>0 && j>0) {
                            float sub_score = substitution_score_prob(start1+i,start2+j);
                            float match_score = mat1(i-1,j-1)   + sub_score + m1_extend;
                            float xgap_score =  mat1_2(i-1,j-1) + sub_score + i1d1_to_m1;
                            float ygap_score =  mat1_3(i-1,j-1) + sub_score + i1d1_to_m1;
                            max_score = maximum(match_score, xgap_score, ygap_score);
                            if(max_score > top_score){
                                top_score = max_score;
                            }
                            if(match_score == max_score) {
                                ts_ptr1(i,j) = match;
                            } else if(xgap_score == max_score) {
                                ts_ptr1(i,j) = xgap;
                            } else if(ygap_score == max_score) {
                                ts_ptr1(i,j) = ygap;
                            } else {
                                cout<<"first no match found, exiting"<<endl;
                                exit(0);
                            }
                            mat1(i,j) = max_score;
                        } else {
                            mat1(i,j) = large_neg;
                            ts_ptr1(i,j) = none;
                        }

                        if(i>0)
                        {
                            float match_open = mat1(i-1,j) + log_emit_indel + m1_to_i1d1;
                            float x_extend = mat1_2(i-1,j) + log_emit_indel + indel_extend;
                            float y_open =   mat1_3(i-1,j) + log_emit_indel + indel_to_indel;
                            max_score = maximum(match_open, x_extend, y_open);
                            if(max_score > top_score){
                                top_score = max_score;
                            }
                            
                            if(match_open == max_score) {
                                ts_ptr1_2(i,j) = match;
                            } else if(x_extend == max_score) {
                                ts_ptr1_2(i,j) = xgap;
                            } else {
                                ts_ptr1_2(i,j) = ygap;
                            }
                            mat1_2(i,j) = max_score;
                        } else {
                            mat1_2(i,j) = large_neg;
                            ts_ptr1_2(i,j) = none;
                        }
                        

                        if(j>0)
                        {   
                            float match_open = mat1(i,j-1) + m1_to_i1d1; // + log_emit_indel;
                            float x_open =   mat1_2(i,j-1) + indel_to_indel; // + log_emit_indel;
                            float y_extend = mat1_3(i,j-1) + indel_extend; // + log_emit_indel;
                            max_score = maximum(match_open, y_extend, x_open);
                            if(max_score > top_score){
                                top_score = max_score;
                            }
                        
                            if(match_open == max_score) {
                                ts_ptr1_3(i,j) = match;
                            } else if(y_extend == max_score) {
                                ts_ptr1_3(i,j) = ygap;
                            } else {
                                ts_ptr1_3(i,j) = xgap;
                            }
                            mat1_3(i,j) = max_score;
                        } else {
                            mat1_3(i,j) = large_neg;
                            ts_ptr1_3(i,j) = none;
                        }
                    }

                    if(local)
                    {
                        mat1(i,j) = 0;
                        mat1_2(i,j) = 0;
                        mat1_3(i,j) = 0;
                        
                        ts_ptr1(i,j) = none;
                        ts_ptr1_2(i,j) = none;
                        ts_ptr1_3(i,j) = none;
                    }
                    if(top_score>maxs)
                    {
                        maxs=top_score;
                        maxi=i;
                        maxj=j;
                    }
                }
            }
            if(debug_matrix) cout<<"\n";
            if(debug_matrix) cout<<"maxs: "<<maxs<<" maxi: "<<maxi<<" maxj: "<<maxj<<"\n";


        for(int i=0;i<=sl1;i++  )
            {

                float pscore1 = large_neg;
                float pscore2 = large_neg;
                float pscore3 = large_neg;
                int pk1 = -1;
                int pk2 = -1;
                int pk3 = -1;

                // from mat1
                for(int k=0;k<=sl2;k++)
                {
                    if(mat1(i-1,k) > pscore1) {
                         pscore1 = mat1(i-1,k);
                         pk1 = k;
                    }
                    if(mat1_2(i-1,k) > pscore2) {
                         pscore2 = mat1_2(i-1,k);
                         pk2 = k;
                    }
                    if(mat1_3(i-1,k) > pscore3) {
                         pscore3 = mat1_3(i-1,k);
                         pk3 = k;
                    }
                }
            // mat2
                for(int j=sl2;j>=0;j--)
                {
                    if(debug_matrix) cout<<start2+j<<" "<<string("-ACGTN").at(seq2.at(start2+j)+1);

                    float score = large_neg;
                    float subst = 0;
                    int ptr = none;
                    int len = 0;
                    float sco = 0;
                    int dist = sl2;
                    int max_jump = 0;
                    float maximum_jump = large_neg;
                    float jump1,jump2,jump3 = large_neg;
                    if(i>0 && j>0)
                    {
                        ptr = match;
                        subst = rev_substitution_score_prob(start1+i,start2+j);
                        
                        if(j<sl2)
                        {
                            score = subst + mat2(i-1,j+1) + m2_extend;
                            len = len2(i-1,j+1)+1;
                            sco = sco2(i-1,j+1) + subst + m2_extend;
                            dist = dist2(i-1,j+1);
                        }
                        jump1 = subst + pscore1 + m1_to_m2;
                        jump2 = subst + pscore2 + i1d1_to_m2;
                        jump3 = subst + pscore3 + i1d1_to_m2;
                        maximum_jump = maximum(jump1,jump2,jump3);
                        // from mat1
                        if(maximum_jump > score)
                        {
                            if(jump1 == maximum_jump) {
                                score = subst + pscore1 + m1_to_m2;
                                ptr = pk1;
                                len = 1;
                                sco = subst;
                                dist = abs(j-pk1);
                                max_jump = 1;
                            } else if(jump2 == maximum_jump) {
                                score = subst + pscore2 + i1d1_to_m2;
                                ptr = pk2;
                                len = 1;
                                sco = subst;
                                dist = abs(j-pk2);
                                max_jump = 2;
                            } else if(jump3 == maximum_jump) {
                                score = subst + pscore3 + i1d1_to_m2;
                                ptr = pk3;
                                len = 1;
                                sco = subst;
                                dist = abs(j-pk3);
                                max_jump = 3;
                            } else {
                                cout<<"none equal max jumpm exiting."<<endl;
                                exit(0);
                            }
                        }
                    }

                    mat2(i,j) = score;
                    ptr2(i,j) = ptr;
                    jump_mat(i,j) = max_jump;
                    len2(i,j) = len;
                    sco2(i,j) = sco;

                    dist2(i,j) = dist;

                }
            }

            maxs = 0;
            maxi = 0;
            maxj = 0;

            for(int i=0;i<=sl1;i++)
            {
                for(int j=0;j<=sl2;j++)
                {
                    if(debug_matrix) cout<<start2+j<<" "<<string("-ACGTN").at(seq2.at(start2+j)+1);
                    // mat3
                    float score = large_neg;
                    int len = 0;
                    int sco = 0;
                    float top_score = large_neg;

                    if(i>0 && j>0)
                    {
                        if(not local && force_overlap && j+start2<clus_start2) {
                            ;
                        } else if(not local && force_overlap && i+start1<clus_start1) {
                            ;
                        } else if(not local && force_overlap && j+start2>true_end2) { 
                            ;
                        }
                        else
                        {
                            float subst = substitution_score_prob(start1+i,start2+j);
                            
                            float match_score =  mat3(i-1, j-1) + subst + m3_extend;
                            float xgap_score =  mat3_2(i-1,j-1) + subst + i3d3_to_m3;
                            float ygap_score = mat3_3(i-1, j-1) + subst + i3d3_to_m3;
                            float max_score = maximum(match_score, xgap_score, ygap_score);
                            if(match_score == max_score) {
                                ptr3(i,j) = match;
                            } else if(xgap_score == max_score) {
                                ptr3(i,j) = xgap;
                            } else if(ygap_score == max_score) {
                                ptr3(i,j) = ygap;
                            } else {
                                cout<<"m3 no equal max, exiting."<<endl;
                                exit(0);
                            }
                            mat3(i,j) = max_score;
                            score = max_score;
                            
                            len = len3(i-1,j-1);
                            

                            float match_open = mat3(i-1,j) + log_emit_indel + m3_to_i3d3;
                            float x_extend = mat3_2(i-1,j) + log_emit_indel + indel_extend;
                            float y_open = mat3_3(i-1,j) +   log_emit_indel + indel_to_indel;
                            float max_score_two = maximum(match_open, x_extend, y_open);

                            if(match_open == max_score_two) {
                                ptr3_2(i,j) = match;
                            } else if(x_extend == max_score_two) {
                                ptr3_2(i,j) = xgap;
                            } else if(y_open == max_score_two) {
                                ptr3_2(i,j) = ygap;
                            } else {
                                cout<<"m2 no equal max2, exiting."<<endl;
                                exit(0);
                            }
                            
                            mat3_2(i,j) = max_score_two;                            

                            match_open =       mat3(i,j-1) + m3_to_i3d3; //+ log_emit_indel;
                            float y_extend = mat3_3(i,j-1) + indel_extend; //+ log_emit_indel;
                            float x_open =   mat3_2(i,j-1) + indel_to_indel; //+ log_emit_indel;
                            float max_score_three = maximum(match_open, y_extend, x_open);

                            if(match_open == max_score_three) {
                                ptr3_3(i,j) = match;
                            } else if(y_extend == max_score_three) {
                                ptr3_3(i,j) = ygap;
                            } else if(x_open == max_score_three) {
                                ptr3_3(i,j) = xgap;
                            } else {
                                cout<<"m3 no equal max3, exiting."<<endl;
                                exit(0);
                            }
                            
                            mat3_3(i,j) = max_score_three;                            
                            
                            top_score = maximum(max_score, max_score_two, max_score_three);
                            
                            if(i==sl1 && j>=true_end2)
                            {
                                max_score_three = maximum(mat3(i,j-1), mat3_2(i,j-1), mat3_3(i,j-1));
                                if(max_score_three > top_score)
                                {
                                    mat3_3(i,j) = max_score_three;
                                    ptr3_3(i,j) = ygap;

                                    len = len3(i,j-1);
                                    sco = sco3(i,j-1);
                                }
                            }

                            // from mat2
                            int dist = sl2;
                            for(int k=0;k<=sl2;k++)
                            {
                                    // jump from m2 to m3
                                    if(subst + m2_to_m3 + mat2(i-1,k) > max_score ||
                                       ( subst + m2_to_m3 + mat2(i-1,k) == max_score && dist2(i-1,k)<dist) )
                                    {
                                        score = subst + mat2(i-1,k) + m2_to_m3;
                                        mat3(i,j) = score;
                                        ptr3(i,j) = k;

                                        len = len2(i-1,k);
                                        sco = sco2(i-1,k);
                                        dist = dist2(i-1,k);

                                        if(score>maxs)
                                        {
                                            maxs=score;
                                            maxi=i;
                                            maxj=k;
                                        }
                                    }
                                    // jump from m2 to i3
                                    if(log_emit_indel + m2_to_i3d3 + mat2(i-1,k) > max_score_two ||
                                       (log_emit_indel + m2_to_i3d3 + mat2(i-1,k) == max_score_two && dist2(i-1,k)<dist) )
                                    {
                                        score = log_emit_indel + mat2(i-1,k) + m2_to_i3d3;
                                        
                                        mat3_2(i,j) = score;
                                        ptr3_2(i,j) = k;

                                        len = len2(i-1,k);
                                        sco = sco2(i-1,k);
                                        dist = dist2(i-1,k);

                                        if(score>maxs)
                                        {
                                            maxs=score;
                                            maxi=i;
                                            maxj=k;
                                        }

                                    }
                                    // jump from m2 to d3
                                    // removed the log_emit_indel from the start of the following two parentheses
                                    // and the score to reflect 0 penalty deletions
                                    if(m2_to_i3d3 + mat2(i-1,k) > max_score_three ||
                                       (m2_to_i3d3 + mat2(i-1,k) == max_score_three && dist2(i-1,k)<dist) )
                                    {
                                        score = mat2(i-1,k) + m2_to_i3d3;
                                        
                                        mat3_3(i,j) = score;
                                        ptr3_3(i,j) = k;

                                        len = len2(i-1,k);
                                        sco = sco2(i-1,k);
                                        dist = dist2(i-1,k);

                                        if(score>maxs)
                                        {
                                            maxs=score;
                                            maxi=i;
                                            maxj=k;
                                        }
                                    } 
                            }
                        }
                    } else if( not local) {
                        mat3(i,j) = large_neg;
                        mat3_2(i,j) = large_neg;
                        mat3_3(i,j) = large_neg;

                        ptr3(i,j) = none;
                        ptr3_2(i,j) = none;
                        ptr3_3(i,j) = none;

                    }

                    if(not local)
                    {
                        len3(i,j) = len;
                        sco3(i,j) = sco;
                    }

                    if(local)
                    {
                        mat3(i,j) = 0;
                        mat3_2(i,j) = 0;
                        mat3_3(i,j) = 0;

                        ptr3(i,j) = none;
                        ptr3_2(i,j) = none;
                        ptr3_3(i,j) = none;

                        len3(i,j) = 0;
                        sco3(i,j) = 0;
                    }
                    
                }
            }

            if(debug_matrix) cout<<maxs<<" "<<maxi<<" "<<maxj<<"\n";
            
            if(debug_score_matrix) {
                cout<<"M_1"<<endl;
                for(int i=0;i<=sl1;i++) {
                    for(int j=0;j<=sl2;j++) {
                        cout<<mat1(i,j)<<",";
                   }
                cout<<endl;
                }
                cout<<"M_2"<<endl;
                for(int i=0;i<=sl1;i++) {
                    for(int j=0;j<=sl2;j++) {
                        cout<<mat2(i,j)<<",";
                    }
                cout<<endl;
                }
                cout<<"M_3"<<endl;
                for(int i=0;i<=sl1;i++) {
                    for(int j=0;j<=sl2;j++) {
                        cout<<mat3(i,j)<<",";
                    }
                cout<<endl;
                }
            exit(0);
            }

            //declare traceback variables
            Array2D<int> *current_ptr = &ptr3;
            int lasti, lastj;
            int current = 0;
            
            //variables to find traceback starting value
            //Array2D<int> *tmp_max_ptr3 = &ptr3;
            Array2D<double> *ts_sco_mat = &mat3;
            Array2D<int> *max_ptr3 = &ptr3;
            int max_match_j = 0;
            float max_match_sco = -10000;
            float max_match_sco_hold = 0;
            int sl1_i = sl1;
            // Determine highest scoring cell at the last position of the descendent
            // sequence, record this cell to get the local alignment score for the
            // template switch alignment, and use this position to initiate traceback.
            
            if(last_match) {
                for(;max_match_j==0;) {
                    for(int j=sl2;j>0;j--) {
                        //current_mat_hold = 0;
                        max_match_sco_hold = maximum(mat3(sl1_i,j), mat3_2(sl1_i,j), mat3_3(sl1_i,j));
                        if(max_match_sco_hold == mat3(sl1_i,j)) {
                            max_ptr3 = &ptr3;
                        } else if(max_match_sco_hold == mat3_2(sl1_i,j)) {
                            max_ptr3 = &ptr3_2;
                        } else if(max_match_sco_hold == mat3_3(sl1_i,j)) {
                            max_ptr3 = &ptr3_3;
                        } else {
                            cout<<"none equal max mat, exiting."<<endl;
                            exit(0);
                        }
                        if((*max_ptr3)(sl1_i,j) == match && max_match_sco_hold > max_match_sco) {
                            max_match_j = j;
                        }
                    }
                }
            } else {
                for(int j=sl2;j>=0;j--) {
                    max_match_sco_hold = maximum(mat3(sl1,j), mat3_2(sl1,j), mat3_3(sl1,j));
                    
                    if(max_match_sco_hold > max_match_sco) {
                        max_match_j = j;
                        if(max_match_sco_hold == mat3(sl1,j)) {
                            current_ptr = &ptr3;
                            ts_sco_mat = &mat3;
                            max_match_sco = mat3(sl1,j);
                        } else if(max_match_sco_hold == mat3_2(sl1,j)) {
                            current_ptr = &ptr3_2;
                            ts_sco_mat = &mat3_2;
                            max_match_sco = mat3_2(sl1,j);
                        } else if(max_match_sco_hold == mat3_3(sl1,j)) {
                            current_ptr = &ptr3_3;
                            ts_sco_mat = &mat3_3;
                            max_match_sco = mat3_3(sl1,j);
                        } else {
                            cout<<"none equal max mat, exiting."<<endl;
                            exit(0);
                        }
                    }
                }
            }
            lasti=sl1;
            lastj=max_match_j;

            traceback_max_sco = maximum(mat3(sl1,sl2), mat3_2(sl1,sl2), mat3_3(sl1,sl2)); 
            ts_sco = (*ts_sco_mat)(lasti,lastj);
            //ts_sco = max_match_sco;
            int i = lasti;
            int j = lastj;
            //current_ptr = max_ptr3;
            current = 1;
            bool ptr_found = false;

            for(;i>=0 || j>=0;)
            {
                seqCoordinate c;
                c.matrix = 3;
                c.pos_x = start1+i;
                c.pos_y = start2+j;
                if(current==1)
                {
                    if((*current_ptr)(i,j) == match) {
                        path->push_back(c);
                        i--;j--;
                        current_ptr = &ptr3;
                        current = 1;
                    } else if((*current_ptr)(i,j) == xgap) {
                        path->push_back(c);
                        i--;j--;
                        current_ptr = &ptr3_2;
                        current = 2;
                    } else if((*current_ptr)(i,j) == ygap) {
                        path->push_back(c);
                        i--;j--;
                        current_ptr = &ptr3_3;
                        current = 3;
                    } else if((*current_ptr)(i,j) == none) {
                        cout<<"i,j: "<<i<<" "<<j<<endl;
                        cout<<"current ptr: "<<(*current_ptr)(i,j)<<endl;
                        cout<<"mat3 i,js: "<<mat3(i,j)<<" "<<mat3_2(i,j)<<" "<<mat3_3(i,j)<<endl;
                        cout<<"ptr3 = none. weird. exiting.\n";
                        exit(0);
                    } else {
                        path->push_back(c);

                        points->at(3).i = start1+i;
                        points->at(3).j = start2+j;
                        

                        j = (*current_ptr)(i,j);
                        i--;

                        points->at(2).i = start1+i;
                        points->at(2).j = start2+j;

                        ptr_found = true;
                        break;
                    }
                    
                }
                
                else if(current==2)
                {
                    if((*current_ptr)(i,j) == match) {
                        c.pos_y = -1;
                        path->push_back(c);
                        i--;
                        current_ptr = &ptr3;
                        current = 1;
                    } else if((*current_ptr)(i,j) == xgap) {
                        c.pos_y = -1;
                        path->push_back(c);
                        i--;
                        current_ptr = &ptr3_2;
                        current = 2;
                    } else if((*current_ptr)(i,j) == ygap) {
                        c.pos_y = -1;
                        path->push_back(c);
                        i--;
                        current_ptr = &ptr3_3;
                        current = 3;
                    } else if((*current_ptr)(i,j) == none) {
                        break;
                    } else {
                        path->push_back(c);
                        points->at(3).i = start1+i;
                        points->at(3).j = start2+j;

                        j = (*current_ptr)(i,j);
                        i--;

                        points->at(2).i = start1+i;
                        points->at(2).j = start2+j;

                        ptr_found = true;
                        break;
                    }
                }

                else if(current==3)
                {
                    if((*current_ptr)(i,j) == match) {
                        c.pos_x = -1;
                        path->push_back(c);
                        j--;
                        current_ptr = &ptr3;
                        current = 1;
                    } else if((*current_ptr)(i,j) == xgap) {
                        c.pos_x = -1;
                        path->push_back(c);
                        j--;
                        current_ptr = &ptr3_2;
                        current = 2;
                    } else if((*current_ptr)(i,j) == ygap) {
                        c.pos_x = -1;
                        path->push_back(c);
                        j--;
                        current_ptr = &ptr3_3;
                        current = 3;
                    } else if((*current_ptr)(i,j) == none) {
                        exit(0);
                    } else {
                        path->push_back(c);

                        points->at(3).i = start1+i;
                        points->at(3).j = start2+j;

                        j = (*current_ptr)(i,j);
                        i--;

                        points->at(2).i = start1+i;
                        points->at(2).j = start2+j;

                        ptr_found = true;
                        break;
                    }
                }

                else
                {
                    cout<<"traceback 1 current not equal 1,2,3"<<endl;
                    cout<<"current: "<<current<<endl;
                    cout<<"i,j: "<<i<<" "<<j<<endl;
                    cout<<(*current_ptr)(i,j)<<endl;
                    exit(0);
                }
            }
            
            if(!ptr_found)
            {
                return;
            }

            ptr_found = false;

            for(;i>=0 || j<=sl2;)
            {
                seqCoordinate c;
                c.matrix = 2;
                c.pos_x = start1+i;
                c.pos_y = start2+j;
                path->push_back(c);

                if(ptr2(i,j)==match)
                {
                    i--;j++;
                }

                else if(ptr2(i,j)==none)
                {
                    cout<<"p2 I J IF NOT FOUND: "<<i<<" "<<j<<" "<<endl;
                    cout<<"ptr2 = none. weird. exiting.\n";
                    exit(0);
                }

                else
                {
                    current = jump_mat(i,j);
                    if(current==1) {
                        current_ptr = &ts_ptr1;
                    } else if(current==2) {
                        current_ptr = &ts_ptr1_2;
                    } else if(current==3) {
                        current_ptr = &ts_ptr1_3;
                    } else {
                        cout<<"current not equal 123, exit."<<endl;
                        cout<<"i,j: "<<i<<" "<<j<<endl;
                        cout<<jump_mat(i,j);
                        exit(0);
                    }

                    points->at(1).i = start1+i;
                    points->at(1).j = start2+j;
                    j = ptr2(i,j);
                    i--;

                    points->at(0).i = start1+i;
                    points->at(0).j = start2+j;

                    ptr_found = true;
                    break;
                }
            }

            if(!ptr_found)
            {
                cout<<"backtracking 2nd path failed. exiting.\n";
                exit(0);
            }
            
            //mat1_traceback
            for(;i>=0 || j>=0;)
            {
                seqCoordinate c;
                c.matrix = 1;
                c.pos_x = start1+i;
                c.pos_y = start2+j;
                
                if(current==1)
                {
                    if((*current_ptr)(i,j) == match) {
                        path->push_back(c);
                        i--;j--;
                        current_ptr = &ts_ptr1;
                        current = 1;
                    } else if((*current_ptr)(i,j) == xgap) {
                        path->push_back(c);
                        i--;j--;
                        current_ptr = &ts_ptr1_2;
                        current = 2;
                    } else if((*current_ptr)(i,j) == ygap) {
                        path->push_back(c);
                        i--;j--;
                        current_ptr = &ts_ptr1_3;
                        current = 3;
                    } else if((*current_ptr)(i,j) == none) {
                        break;
                    } else {
                        cout<<"ptr 1 traceback not correct"<<endl;
                        cout<<"i,j: "<<i<<" "<<j<<endl;
                        cout<<"ptr val: "<<(*current_ptr)(i,j)<<endl;
                        exit(0);
                    }

                }
                
                else if(current==2)
                {
                    if((*current_ptr)(i,j) == match) {
                        c.pos_y = -1;
                        path->push_back(c);
                        i--;
                        current_ptr = &ts_ptr1;
                        current = 1;
                    } else if((*current_ptr)(i,j) == xgap) {
                        c.pos_y = -1;
                        path->push_back(c);
                        i--;
                        current_ptr = &ts_ptr1_2;
                        current = 2;
                    } else if((*current_ptr)(i,j) == ygap) {
                        c.pos_y = -1;
                        path->push_back(c);
                        i--;
                        current_ptr = &ts_ptr1_3;
                        current = 3;
                    } else if((*current_ptr)(i,j) == none) {
                        break;
                    } else {
                        cout<<"ptr 1 traceback not correct"<<endl;
                        cout<<"i,j: "<<i<<" "<<j<<endl;
                        cout<<"ptr val: "<<(*current_ptr)(i,j)<<endl;
                        exit(0);
                    }
                }

                else if(current==3)
                {
                    if((*current_ptr)(i,j) == match) {
                        c.pos_x = -1;
                        path->push_back(c);
                        j--;
                        current_ptr = &ts_ptr1;
                        current = 1;
                    } else if((*current_ptr)(i,j) == xgap) {
                        c.pos_x = -1;
                        path->push_back(c);
                        j--;
                        current_ptr = &ts_ptr1_2;
                        current = 2;
                    } else if((*current_ptr)(i,j) == ygap) {
                        c.pos_x = -1;
                        path->push_back(c);
                        j--;
                        current_ptr = &ts_ptr1_3;
                        current = 3;
                    } else if((*current_ptr)(i,j) == none) {
                        break;
                    } else {
                        cout<<"ptr 1 traceback not correct"<<endl;
                        cout<<"i,j: "<<i<<" "<<j<<endl;
                        cout<<"ptr val: "<<(*current_ptr)(i,j)<<endl;
                        exit(0);
                    }
                }
                
                else
                {
                    cout<<"current mat1 traceback not 1,2,3"<<endl;
                    cout<<"current: "<<current<<endl;
                    exit(0);
                }
            }
        }

        /********************    alignment itself   ******************************/



        /********************     alignment scan    ******************************/

        void scan_alignment(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path)
        {
            int width = arg_get("scan-window-width").as<int>();
            int limit = arg_get("scan-window-limit").as<int>();
            int flank = arg_get("scan-flank").as<int>();

            int sum = 0;
            int ss = 0;
            if(arg_is("scan-start"))
            {
                int s = arg_get("scan-start").as<int>();
                if(s>0 && s<slg)
                    ss=s;
            }
            int se = slg;
            if(arg_is("scan-end"))
            {
                int s = arg_get("scan-end").as<int>();
                if(s>ss && s<slg)
                    se=s;
            }

            if(not verbose)
                cout<<"chrom,clus_start_chrom,clus_start_align,clust_start1,clust_end1,sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref,iden_up,ident_rep,ident_down,ident_inv,ident_fwd,ident_epo,masked,sum_ins,sum_del,sum_mis,sum_nuc,CpG,clus_ins,clus_del,clus_mis,fwd_score,ts_score_local,ts_score_global,ts_ref_seq_len,ts_qry_seq_len,frag_L1_size,frag_23_size,frag_4R_size\n";

            int max_length = arg_get("max-event-length").as<int>();
            vector<pair<int,int> > clusters;


            int i = ss;

            int c_first=-1, c_start=-1, c_end=-1;

            // width is the scan window e.g 10bp
            for(;i<width;i++)
            {
                if(fseq1.at(i) != fseq2.at(i))
                {
                    sum++;
                    // set first position of cluster
                    if(sum==1)
                    {
                        c_first = i;
                    }
                }

                 
                if(sum==limit)
                {
                    c_start = i-limit+2;
                    if(c_first < c_start)
                        c_start = c_first;

                }
            }

            for(;i<se;i++)
            {
                // decrease cluster mutatation count if alignment position i-window size
                // is not a match
                if(fseq1.at(i-width) != fseq2.at(i-width) && sum>0)
                    sum--;
                if(fseq1.at(i) != fseq2.at(i))
                {
                    sum++;
                    if(sum==1)
                    {
                        c_first = i;
                    }
                }
                // if minimum number of diffs found and i-window size is a match, define cluster
                if(sum>=limit)
                {
                    int start_pos1 = index1.at(0);
                    int start_pos2 = index2.at(0);

                    // if the first position in the cluster is not at 0, assign cluster start
                    if(i-flank-limit>0)
                    {
                        start_pos1 = index1.at(i-flank-limit);
                        start_pos2 = index2.at(i-flank-limit);
                    }

                    int end_pos1 = index1.at(i);
                    int end_pos2 = index2.at(i);


                    c_start = i-limit+1;
                    if(c_first < c_start)
                        c_start = c_first;

                    clus_start1 = index1.at(c_start);
                    clus_start2 = index2.at(c_start);

                    if(c_start>0)
                    {
                        clus_start1 = index1.at(c_start-1);
                        clus_start2 = index2.at(c_start-1);
                    }

                    i++;
                    // expand cluster outwards from the defined cluster start position
                    for(;i<se;i++)
                    {
                        // decrease cluster mutatation count if alignment position i-window size
                        // is not a match
                        if(fseq1.at(i-width) != fseq2.at(i-width) && sum>0)
                            sum--;
                        if(fseq1.at(i) != fseq2.at(i))
                            sum++;
                        
                        // set the end cluster position to the final position of the input if reached
                        if(i-width+flank<slg)
                            end_pos1 = index1.at(i-width+flank);
                        else
                            end_pos1 = index1.at(slg-1);

                        if(i-width+flank<slg)
                           end_pos2 = index2.at(i-width+flank);
                        else
                            end_pos2 = index2.at(slg-1);

                        // define cluster end if end of input sequence reached 
                        if(sum==0)
                        {

                            c_end = i-width+1;

                            clus_end1 = index1.at(c_end);
                            clus_end2 = index2.at(c_end);

                            break;
                        }
                    }
                    // define cluster end 
                    if(sum>0)
                    {

                        c_end = i-width+1;

                        clus_end1 = index1.at(c_end);
                        clus_end2 = index2.at(c_end);
                    }

                    if(clus_end1-clus_start1<=max_length && clus_end2-clus_start2<=max_length && clus_start1<=clus_end1)
                    {
                        clusters.push_back(make_pair(c_start,c_end));
                        this->check_scan_position(path,points,fwd_path,start_pos1,end_pos1,start_pos2,end_pos2);
                    }
                }
            }

            if(arg_is("cluster-annotation"))
            {
                stringstream str1;
                stringstream str2;
                stringstream anno;
                stringstream repe;

                clusters.push_back(make_pair(se,se));
                vector<pair<int,int> >::iterator it = clusters.begin();

                string alpha = "-ACGTN";
                int i=ss;
                for(;it!=clusters.end();it++)
                {

                    for(;i<it->first;i++)
                    {
                        str1<<alpha.at(fseq1.at(i)+1);
                        str2<<alpha.at(fseq2.at(i)+1);
                        anno<<"A";
                        if(i<(int)index1.size() && index1.at(i)<(int)mask1.size())
                            repe<<alpha.at(mask1.at(index1.at(i))+3);
                        else
                            repe<<"G";
                    }

                    i = it->first;
                    for(;i<it->second;i++)
                    {
                        str1<<alpha.at(fseq1.at(i)+1);
                        str2<<alpha.at(fseq2.at(i)+1);
                        anno<<"C";
                        if(i<(int)index1.size() && index1.at(i)<(int)mask1.size())
                            repe<<alpha.at(mask1.at(index1.at(i))+3);
                        else
                            repe<<"G";
                    }

                    i = it->second;
                }
                for(;i<se;i++)
                {
                    str1<<alpha.at(fseq1.at(i)+1);
                    str2<<alpha.at(fseq2.at(i)+1);
                    anno<<"A";
                    if(i<(int)index1.size() && index1.at(i)<(int)mask1.size())
                        repe<<alpha.at(mask1.at(index1.at(i))+3);
                    else
                        repe<<"G";
                }
                ofstream fout(arg_get("cluster-annotation").as<string>().c_str());
                fout<<">"<<qry_name<<endl<<str1.str()<<endl<<">"<<ref_name<<endl<<str2.str()<<endl<<">cluster"<<endl<<anno.str()<<endl<<">repeat"<<endl<<repe.str()<<endl;
            }
        }



    void check_scan_position(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path,int min1,int max1,int min2,int max2)
    {
        int max_length = arg_get("max-event-length").as<int>();

        start1 = min1;
        start2 = min2;
        end1 = max1;
        end2 = max2;
        true_start2 = min2;
        true_end2 = max2;
        if(end1-start1>max_length || end2-start2>max_length)
            return;

        int flank = arg_get("ref-flank").as<int>();
        // +/100 to ref flank
        if(flank>0 && start2-flank>=0)
            start2 = start2-flank;
        if(flank>0 && end2+flank<(int)seq2.size())
            end2 = end2+flank;

        // length of alignment fragments undergoing re-alignment
        sl1 = end1-start1;
        sl2 = end2-start2;
       
        path->clear();
        fwd_path->clear();
        this->align_sequences(path,points,false);
        this->fwd_align_sequences(fwd_path,path);
        if(!arg_is("long-output"))
        {
            this->print_switch_process(path,points);
            this->print_inversion_fragment(path,points,fwd_path);
        }
        else
        {
            this->print_switch_process(path,points);
            this->print_inversion_fragment(path,points,fwd_path);
        }
    }

    /********************     alignment scan    ******************************/




    /********************     print events      ******************************/

    void print_events(vector<seqCoordinate> *path,vector<switchPoint> *points,vector<seqCoordinate> *fwd_path)
    {
        verbose = true;
        int flank = arg_get("scan-flank").as<int>();
        string filename = arg_get("print-file").as<string>();

        ifstream input(filename.c_str(), ios::in);
        if (!input) { return; }

        string temp;
        getline(input, temp, '\n');
        
        // Previous line both caused errors and didn't print X/Y chroms
        // without an alphanum line appended before, fixed.
        bool has_line = temp.empty();

        while(!input.eof())
        {
            if(has_line)
                getline(input, temp, '\n');
            size_t field_count = count(temp.begin(), temp.end(), ',');
            
            // If there aren't 32 fields in the csv line, break
            // means only clean lines can be included, this var will need changing if
            // additional output fields added
            if(field_count!=32)
            {
                break;
            }

            has_line = false;


            string chrom = temp.substr(0,temp.find(','));
            temp = temp.substr(temp.find(',')+1);

            if(chrom=="scan finished")
                break;

            stringstream tempstr(temp);

            int chrom_site,aligned,clus_st1,clus_en1,sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref;
            char sep;
            tempstr>>chrom_site>>sep>>aligned>>sep>>clus_st1>>sep>>clus_en1>>sep>>sp1_qry>>sep>>sp1_ref>>sep>>sp2_ref>>sep>>sp3_ref>>sep>>sp4_ref>>sep>>temp;

            // comment out to clean output
            //if(arg_is("long-output"))
            //    cout<<chrom<<","<<tempstr.str()<<endl;

            // Procedure below adds e.g. 40nt to start/end of cluster coords
            // in each sequence. This will yield different length sequences as
            // the cluster is identified as the same size in the alignment but
            // the e.g. +/-40nt alignment region surrounding this may contain
            // indels.
            int start_pos1 = index1.at(0);
            int start_pos2 = index2.at(0);
            if(clus_st1>flank)
            {
                start_pos1 = index1.at(rindex1.at(clus_st1)-flank);
                start_pos2 = index2.at(rindex1.at(clus_st1)-flank);
            }

            int end_pos1 = index1.at(rindex1.at(clus_en1));
            int end_pos2 = index2.at(rindex1.at(clus_en1));
            if(rindex1.at(clus_en1)+flank<(int)index1.size())
                end_pos1 = index1.at(rindex1.at(clus_en1)+flank);
            else
                end_pos1 = index1.at(index1.size()-1);

            if(rindex1.at(clus_en1)+flank<(int)index2.size())
                end_pos2 = index2.at(rindex1.at(clus_en1)+flank);
            else
                end_pos2 = index2.at(index2.size()-1);

            clus_start1 = clus_st1;
            clus_start2 = index2.at(rindex1.at(clus_st1));

            clus_end1 = clus_en1;
            clus_end2 = index2.at(rindex1.at(clus_en1));
            
            this->check_scan_position(path,points,fwd_path,start_pos1,end_pos1,start_pos2,end_pos2);

        }
    }


   /********************     print events      ******************************/


public:
    TSApHMM() {}

    int run_alignment(int argc, char *argv[])
    {

        Fasta_entry ref;
        Fasta_entry qry;

        pair_data = false;
        maximise_score = false;
        maximise_length = false;
        maximise_length_alt = false;
        verbose = false;
        
        string ref_chrom;
        ref_chrom = "";
        int ref_chrom_len;
        qry_name = "";
        ref_name = "";

        try {
            this->read_command_line_arguments(argc,argv);

            chrom = "0";
            chrom_start = 0;

            if(arg_is("pair"))
            {
                string file = arg_get("pair").as<string>();

                this->read_pair_data(&file,&qry,&ref);

                string temp = qry.name;
                // get value in reference chromosome option
                if(arg_is("ref_chrom"))
                    ref_chrom = arg_get("ref_chrom").as<string>();
                ref_chrom_len = ref_chrom.length();
                if(temp.substr(0,ref_chrom_len)==ref_chrom)
                {
                    size_t pos1 = temp.find_first_of('.');
                    size_t pos2 = temp.find_first_of('.',pos1+1);
                    chrom = temp.substr(pos1+1,pos2-pos1-1);
                    pos1=pos2;
                    pos2 = temp.find_first_of('.',pos1+1);
                    chrom_start = atoi(temp.substr(pos1+1,pos2-pos1-1).c_str());
                }
                else
                {
                    temp = ref.name;
                    if(temp.substr(0,ref_chrom_len)==ref_chrom)
                    {
                        size_t pos1 = temp.find_first_of('.');
                        size_t pos2 = temp.find_first_of('.',pos1+1);
                        chrom = temp.substr(pos1+1,pos2-pos1-1);
                        pos1=pos2;
                        pos2 = temp.find_first_of('.',pos1+1);
                        chrom_start = atoi(temp.substr(pos1+1,pos2-pos1-1).c_str());
                    }
                }
            }
            else
            {
                if(arg_is("ref"))
                {
                    string file = arg_get("ref").as<string>();
                    this->read_data(&file,&ref);
                }

                if(arg_is("qry"))
                {
                    string file = arg_get("qry").as<string>();
                    this->read_data(&file,&qry);
                }
            }

            if(ref.length==0 || qry.length==0)
            {
                cout<<"Error reading the sequence input. Exiting."<<endl;
                return 1;
            }
            if(arg_is("qry-name"))
                qry_name = arg_get("qry-name").as<string>();
            if(arg_is("ref-name"))
                ref_name = arg_get("ref-name").as<string>();


        } catch ( const po::error& e ) {
            this->arg_error();
        }

        vector<seqCoordinate> path;
        vector<switchPoint> points;

        for(int i=0;i<4;i++)
        {
            switchPoint sp;
            points.push_back(sp);
        }

        vector<seqCoordinate> fwd_path;

        if(verbose)
        {
            cout<<"\nQuery:     "<<qry.name<<endl;
            cout<<"Reference: "<<ref.name<<endl<<endl;
        }

        if(arg_is("pair"))
        {
            this->build_indeces(&qry.sequence,&ref.sequence);
            this->build_sequences(&qry.sequence,&ref.sequence);

            if(arg_is("scan"))
            {
                this->scan_alignment(&path,&points,&fwd_path);
                cout<<"scan finished\n";
            }
            else if(arg_is("print-file"))
            {
                this->print_events(&path,&points,&fwd_path);
            }
            else
            {
                this->set_alignment_fragments();
                this->align_sequences(&path,&points);
                this->fwd_align_sequences(&fwd_path,&path);
                this->print_inversion_fragment(&path,&points,&fwd_path);
                this->print_switch_process(&path,&points);
            }
        }
        else
        {
            this->build_sequences(&qry.sequence,&ref.sequence);
            this->set_two_fragments();

            this->align_sequences(&path,&points);
            this->fwd_align_sequences(&fwd_path,&path);
            this->print_inversion_fragment(&path,&points,&fwd_path);
            this->print_switch_process(&path,&points);
        }

        return 0;
    }
};

int main(int argc, char *argv[])
{
    TSApHMM tsa_pairhmm;
    return tsa_pairhmm.run_alignment(argc,argv);
}
