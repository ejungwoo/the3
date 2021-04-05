#ifndef __summary2_binning_hh__
#define __summary2_binning_hh__

#include "TH1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TAxis.h"
#include <iostream>
using namespace std;

class binning;
class binning2;

class binning
{
  public:
    int fN = 0; ///< number of bins
    double fMin = 0; ///< lower bound
    double fMax = 0; ///< upper bound
    double fW = 0; ///< binning space width
    double fValue = 0; ///< value will be set after iteration using next(), back(), nextb()
    int fIdx = 0; ///< index for iteration
    const char *fTitle = "";
    const char *fExpression = "";

  public:
    binning(int n, double min, double max, const char *ttl="", const char *expr="");
    ~binning() {}

    void init();
    void make(TTree *tree, const char* branchName);

    bool   isNull()               const { if (fN<1||(fMin==0&&fMax==0)) return true; return false; }
    int    getN()                 const { return fN; }
    double getMin()               const { return fMin; }
    double getMax()               const { return fMax; }
    double getW()                 const { return fW; }
    double getValue()             const { return fValue; }
    const char *getTitle()        const { return fTitle; }
    const char *getExpression()   const { return fExpression; }
    const char *getE()            const { return fExpression; }

    int    getIdx()               const { return fIdx; }
    int    getBin (double val)    const { return int((val-fMin)/fW); }
    int    findBin(double val)    const { return getBin(val); }

    //iterator
    void   reset(bool iOU=0);
    bool   next (bool iOU=0);

    void   end();
    bool   prev();

    int    ii()     const { return fIdx-1; }
    int    bi()     const { return fIdx; }
    double cValue()   const { return fValue; }
    double cLow()     const { return lowEdge(fIdx); }
    double cHigh()    const { return highEdge(fIdx); }

    double getFullWidth()         const { return .5*(fMax + fMin); }
    double getFullCenter()        const { return (fMax - fMin); }
    double lowEdge  (int bin= 1)  const { return fMin+(bin-1)*(fMax-fMin)/fN; }
    double highEdge (int bin=-1)  const { if (bin==-1) bin=fN; return fMin+(bin)*(fMax-fMin)/fN; }
    double getCenter(int bin)     const { return (fMin + (bin-.5)*fW); }
    bool   isInside(double value) const { if(value>fMin && value<fMax) return true; return false; }

    void set     (double n, double min, double max, const char *ttl="", const char *expression="");
    void setN    (double n)          { fN = n; fW = (fMax-fMin)/fN; }
    void setW    (double w)          { fW = w; fN = int((fMax-fMin)/fW); }
    void setMin  (double min)        { fMin = min; fW = (fMax-fMin)/fN; }
    void setMax  (double max)        { fMax = max; fW = (fMax-fMin)/fN; }
    void setTitle(const char *ttl)   { fTitle = ttl; }
    void setExpression(const char *expression)   { fExpression = expression; }

    TH1D *newHist(const char *name, const char *title="");

    TString print(bool pout=1) const;
    void operator=(const binning binn);

    binning2 operator*(const binning binn);
};


class binning2 {
  public :
    const char *fTitleX;
    int fNX = 0;
    double fMinX = 0;
    double fMaxX = 0;
    const char *fTitleY;
    int fNY = 0;
    double fMinY = 0;
    double fMaxY = 0;

    binning2(const char *titleX,int nX, double minX, double maxX, const char *titleY, int nY, double minY, double maxY)
    : fTitleX(titleX) ,fNX(nX), fMinX(minX), fMaxX(maxX), fTitleY(titleY) ,fNY(nY), fMinY(minY), fMaxY(maxY) {}

    TH2D *newHist(const char *nameHist, const char *titleHist="") {
      auto titlexy = Form("%s;%s;%s;",titleHist,fTitleX,fTitleY);
      auto hist = new TH2D(nameHist,titlexy,fNX,fMinX,fMaxX,fNY,fMinY,fMaxY);
      return hist;
    }
};



binning::binning(int n, double min, double max, const char *ttl, const char *expression) : fN(n), fMin(min), fMax(max), fTitle(ttl), fExpression(expression) { init(); }

void binning::reset(bool includeOUFlow)  { fIdx = (includeOUFlow?-1:0); }

bool binning::next(bool includeOUFlow) 
{
  if ((includeOUFlow && fIdx>fN) || (fIdx>fN-1))
    return false;
  fValue = fMin + (fIdx++) * fW + .5 * fW;
  return true;
}

void binning::end() { fIdx = fN+1; }

bool binning::prev()
{
  if (fIdx<2)
    return false;

  fValue = fMin + (fIdx--) * fW + .5 * fW;
  return true;
}


TH1D *binning::newHist(const char *name, const char *title) {
  auto titlexy = Form("%s;%s;",title,fTitle);
  return (new TH1D(name,titlexy,fN,fMin,fMax));
}

void binning::init()
{
       if (fW>0&&fN<1) setW(fW);
  else if (fN>0&&fW<1) setN(fN);
}


void binning::set(double n, double min, double max, const char *ttl, const char *expression)
{
  fN = n;
  fMin = min;
  fMax = max;
  fW = (fMax-fMin)/fN;
  if (TString(ttl).IsNull()==false) fTitle = ttl;
  if (TString(expression).IsNull()==false) fExpression = expression;
}

TString binning::print(bool pout) const
{
  auto rm0 = [](TString vstring) {
    auto posdot = vstring.Index(".");
    auto poslast = vstring.Sizeof()-2;
    auto lstring = TString(vstring(poslast));
    if (posdot>=0&&posdot<poslast)
      while (lstring=="0") {
        vstring.Remove(poslast);
        poslast = vstring.Sizeof()-2;
        lstring = TString(vstring(poslast));
      }
    return vstring;
  };

  TString minString = rm0(Form("%f",fMin));
  TString maxString = rm0(Form("%f",fMax));
  TString wString = rm0(Form("%f",fW));

  TString line = Form("%d,%s,%s",fN,minString.Data(),maxString.Data());
  if (pout)
    cout << line << endl;
  return line;
}

void binning::operator=(const binning binn)
{
  fN = binn.getN();
  fMin = binn.getMin();
  fMax = binn.getMax();
  fW = binn.getW();
  fTitle = binn.getTitle();
}

binning2 binning::operator*(const binning binn) {
  return binning2(fTitle,fN,fMin,fMax, binn.getTitle(),binn.getN(),binn.getMin(),binn.getMax());
}

void binning::make(TTree *tree, const char* branchName)
{
  if (fN<1) fN = 100;
  if (tree==nullptr) {
    fMin = 0;
    fMax = 1;
    return;
  }
  auto min = tree -> GetMinimum(branchName);
  auto max = tree -> GetMaximum(branchName);
  if (0) {
    auto dmm = (max - min) / 10.;
    fMax = max + dmm;
    fMin = min - dmm;
  }
  else {
    fMax = max;
    fMin = min;
  }
  fW = (fMax-fMin)/fN;
}

#endif
