#include "bwa_fmd_index.hpp"
#include "bwt.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr uint8_t A=0,C=1,G=2,T=3;

uint8_t code(char c) {
    switch (c) {
    case 'A': return A; case 'C': return C; case 'G': return G; case 'T': return T;
    default: throw std::runtime_error("Invalid DNA base: " + std::string(1,c));
    }
}
char base(uint8_t c) {
    static constexpr char b[]={'A','C','G','T'};
    if(c>3) throw std::runtime_error("Invalid base code");
    return b[c];
}
uint8_t comp(uint8_t c) { return static_cast<uint8_t>(3-c); }
char comp(char c) { return base(comp(code(c))); }

std::string rc(const std::string& s) {
    std::string r; r.reserve(s.size());
    for(auto i=s.rbegin(); i!=s.rend(); ++i) r.push_back(comp(*i));
    return r;
}

std::string read_fasta(const std::string& fn) {
    std::ifstream in(fn);
    if(!in) throw std::runtime_error("Cannot open FASTA: "+fn);
    std::string s,line;
    while(std::getline(in,line)) {
        if(line.empty() || line[0]=='>') continue;
        for(char c:line) {
            if(c==' '||c=='\r'||c=='\t') continue;
            if(c>='a'&&c<='z') c=char(c-'a'+'A');
            code(c); s.push_back(c);
        }
    }
    if(s.empty()) throw std::runtime_error("FASTA contains no sequence");
    return s;
}

std::vector<uint64_t> occurrences(const std::string& text,
                                   const std::string& p) {
    std::vector<uint64_t> v;
    if(p.empty() || p.size()>text.size()) return v;
    for(size_t i=0;i+p.size()<=text.size();++i)
        if(text.compare(i,p.size(),p)==0) v.push_back(i);
    return v;
}

std::string pos_string(const std::vector<uint64_t>& v) {
    std::string s="[";
    for(size_t i=0;i<v.size();++i) {
        if(i) s+=", ";
        s+=std::to_string(v[i]);
    }
    return s+"]";
}
std::string iv(const Interval& x) {
    return "["+std::to_string(x.l)+","+std::to_string(x.r)+")";
}
[[noreturn]] void fail(const std::string& s) {
    std::cerr<<"\nFAIL: "<<s<<"\n";
    std::abort();
}
bool same(const SA_Range& a,const SA_Range& b) {
    auto ap=a.primary_interval(), aq=a.companion_interval();
    auto bp=b.primary_interval(), bq=b.companion_interval();
    return ap.l==bp.l&&ap.r==bp.r&&aq.l==bq.l&&aq.r==bq.r;
}

/*
 * The BWA FMD index is built on:
 *
 *     F = R + RC(R)
 *
 * Consequently primary(P) contains occurrences of P in F, while
 * companion(P) contains occurrences of RC(P) in F.
 */
void check_range(const BwaFMDIndex& idx, const SA_Range& r,
                 const std::string& p,
                 const std::vector<uint64_t>& ep,
                 const std::vector<uint64_t>& ec,
                 const std::string& where) {
    auto pi=r.primary_interval(), qi=r.companion_interval();
    if(pi.size()!=ep.size())
        fail(where+" primary size, pattern="+p+
             " actual="+std::to_string(pi.size())+
             " expected="+std::to_string(ep.size()));
    if(qi.size()!=ec.size())
        fail(where+" companion size, pattern="+p+
             " actual="+std::to_string(qi.size())+
             " expected="+std::to_string(ec.size()));
    if(pi.size()!=qi.size())
        fail(where+" primary/companion sizes differ, pattern="+p);

    if(!ep.empty()) {
        auto a=idx.locate(r), e=ep;
        std::sort(a.begin(),a.end());
        std::sort(e.begin(),e.end());
        if(a!=e)
            fail(where+" SA coordinates, pattern="+p+
                 "\n expected="+pos_string(e)+
                 "\n actual  ="+pos_string(a));
    }
}

SA_Range build_right(const BwaFMDIndex& idx,const std::string& p) {
    SA_Range r=idx.initial_range(code(p[0]));
    for(size_t i=1;i<p.size();++i) r=idx.extend_right(r,code(p[i]));
    return r;
}

SA_Range build_left(const BwaFMDIndex& idx,const std::string& p) {
    size_t n=p.size();
    SA_Range r=idx.initial_range(code(p[n-1]));
    for(size_t i=n-1;i-- > 0;) r=idx.extend_left(r,code(p[i]));
    return r;
}
SA_Range build_mixed(const BwaFMDIndex& idx,const std::string& p) {
    size_t m=p.size()/2;
    SA_Range r=idx.initial_range(code(p[m]));
    for(size_t i=m;i-- > 0;) r=idx.extend_left(r,code(p[i]));
    for(size_t i=m+1;i<p.size();++i) r=idx.extend_right(r,code(p[i]));
    return r;
}

void check_rc_symmetry(const BwaFMDIndex& idx,const std::string& F,
                       const std::string& p,const SA_Range& r) {
    const std::string q=rc(p);
    const SA_Range rr=build_right(idx,q);
    auto ep=occurrences(F,p), er=occurrences(F,q);
    check_range(idx,r,p,ep,er,"search(P)");
    check_range(idx,rr,q,er,ep,"search(RC(P))");

    auto a=r.primary_interval(), b=r.companion_interval();
    auto c=rr.primary_interval(), d=rr.companion_interval();
    if(a.l!=d.l||a.r!=d.r||b.l!=c.l||b.r!=c.r)
        fail("RC interval symmetry, pattern="+p+
             "\n primary(P)="+iv(a)+" companion(P)="+iv(b)+
             "\n primary(RC)="+iv(c)+" companion(RC)="+iv(d));
}

void check_all(const BwaFMDIndex& idx,const std::string& F,
               const std::string& p,const SA_Range& r) {
    SA_Range L[4],R[4];
    idx.extend_left_all(r,L);
    idx.extend_right_all(r,R);

    for(uint8_t c=0;c<4;++c) {
        std::string lp=std::string(1,base(c))+p;
        std::string rp=p+std::string(1,base(c));
        auto el=occurrences(F,lp), elc=occurrences(F,rc(lp));
        auto er=occurrences(F,rp), erc=occurrences(F,rc(rp));
        check_range(idx,L[c],lp,el,elc,"extend_left_all");
        check_range(idx,R[c],rp,er,erc,"extend_right_all");
    }
}

void check_singleton(const BwaFMDIndex& idx,const std::string& F,
                     const std::string& p,const SA_Range& r) {
    if(r.size()!=1) return;

    for(uint8_t c=0;c<4;++c) {
        std::string lp=std::string(1,base(c))+p;
        std::string rp=p+std::string(1,base(c));
        auto el=occurrences(F,lp), er=occurrences(F,rp);
        SA_Range fl,fr;
        bool lok=idx.extend_left_singleton(r,c,fl);
        bool rok=idx.extend_right_singleton(r,c,fr);

        if(lok!=(el.size()==1))
            fail("left singleton availability, "+lp);
        if(rok!=(er.size()==1))
            fail("right singleton availability, "+rp);

        if(lok) {
            auto g=idx.extend_left(r,c);
            if(!same(fl,g)) fail("left singleton != general, "+lp);
            check_range(idx,fl,lp,el,occurrences(F,rc(lp)),"left singleton");
        }
        if(rok) {
            auto g=idx.extend_right(r,c);
            if(!same(fr,g)) fail("right singleton != general, "+rp);
            check_range(idx,fr,rp,er,occurrences(F,rc(rp)),"right singleton");
        }
    }
}

void check_halves(const BwaFMDIndex& idx,const std::string& R,
                  const std::string& F,const std::string& p) {
    auto r=build_right(idx,p);
    auto actual=idx.locate(r);
    auto expected_f=occurrences(R,p);
    auto expected_rc=occurrences(rc(R),p);
    std::set<uint64_t> af,ar,ef(expected_f.begin(),expected_f.end()),er;
    for (auto x : actual) {
        if (x + p.size() <= R.size()) {
            af.insert(x);
        } else if (x >= R.size() &&
                   x + p.size() <= F.size()) {
            ar.insert(x);
        }
    }
    for(auto x:expected_rc) er.insert(static_cast<uint64_t>(R.size())+x);
    if(af!=ef) fail("forward-half coordinates, "+p);
    if(ar!=er) fail("RC-half coordinates, "+p);
    (void)F;
}

template<class Fn>
void enumerate(unsigned n,std::string& p,Fn&& fn) {
    if(p.size()==n) { fn(p); return; }
    for(uint8_t c=0;c<4;++c) {
        p.push_back(base(c));
        enumerate(n,p,fn);
        p.pop_back();
    }
}

void test_pattern(const BwaFMDIndex& idx,const std::string& R,
                  const std::string& F,const std::string& p) {
    auto ep=occurrences(F,p);
    auto ec=occurrences(F,rc(p));

    auto a=build_right(idx,p);
    auto b=build_left(idx,p);
    auto c=build_mixed(idx,p);

    check_range(idx,a,p,ep,ec,"build_right");
    check_range(idx,b,p,ep,ec,"build_left");
    check_range(idx,c,p,ep,ec,"build_mixed");

    if(!same(a,b)) fail("left/right construction differs, "+p);
    if(!same(a,c)) fail("mixed construction differs, "+p);

    check_rc_symmetry(idx,F,p,a);
    check_all(idx,F,p,a);
    check_singleton(idx,F,p,a);

    if(!ep.empty())
        check_halves(idx,R,F,p);
}

} // namespace

int main(int argc,char** argv) {
    if(argc!=3) {
        std::cerr<<"Usage: "<<argv[0]<<" <bwa-index-prefix> <reference.fa>\n";
        return 2;
    }

    const std::string prefix=argv[1];
    const std::string bwt_filename=prefix+".bwt";
    const std::string sa_filename=prefix+".sa";
    const std::string fasta_filename=argv[2];

    const std::string R=read_fasta(fasta_filename);
    const std::string RCr=rc(R);
    const std::string F=R+RCr;

    std::cout<<"Reference length: "<<R.size()<<"\n"
             <<"Reference:       "<<R<<"\n"
             <<"Reverse comp:    "<<RCr<<"\n"
             <<"Indexed text:    "<<F<<"\n"
             <<"BWT:             "<<bwt_filename<<"\n"
             <<"SA:              "<<sa_filename<<"\n\n";

    bwt_t* bwt=bwt_restore_bwt(bwt_filename.c_str());
    if(!bwt) {
        std::cerr<<"Could not load BWT: "<<bwt_filename<<"\n";
        return 1;
    }

    /* locate() needs the actual .sa file, not the .bwt filename. */
    bwt_restore_sa(sa_filename.c_str(),bwt);
    BwaFMDIndex idx(bwt);

    try {
        std::cout<<"Checking initial states...\n";
        for(uint8_t c=0;c<4;++c) {
            std::string p(1,base(c));
            check_range(idx,idx.initial_range(c),p,
                        occurrences(F,p),occurrences(F,rc(p)),
                        "initial_range");
        }
        std::cout<<"Initial states passed.\n\n";

        constexpr unsigned MAX_LENGTH=8;
        uint64_t count=0;

        for(unsigned n=1;n<=MAX_LENGTH;++n) {
            std::string p;
            p.reserve(n);
            enumerate(n,p,[&](const std::string& x) {
                test_pattern(idx,R,F,x);
                ++count;
            });
            std::cout<<"Length "<<n<<": passed\n";
        }

        std::cout<<"\n========================================\n"
                 <<"All exhaustive BWA FMD tests passed.\n"
                 <<"Patterns tested: "<<count<<"\n"
                 <<"Indexed text: R + RC(R)\n"
                 <<"========================================\n";
    } catch(...) {
        bwt_destroy(bwt);
        throw;
    }

    bwt_destroy(bwt);
    return 0;
}
