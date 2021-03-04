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

  public:
    binning(int n=-1, double min=0, double max=0, double w=-1);
    binning(int n, double min, double max, const char *ttl);
    binning(binning const & binn) : fN(binn.fN), fMin(binn.fMin), fMax(binn.fMax), fW((fMax-fMin)/fN) {}
    binning(TH1 *hist, int i=1);
    binning(TGraph *graph, int i=1);
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

    int    getIdx()               const { return fIdx; }
    int    getBin (double val)    const { return int((val-fMin)/fW); }
    int    findBin(double val)    const { return getBin(val); }

    //iterator
    void   reset(bool iOU=0);
    bool   next (bool iOU=0);
    int    cIdx()     const { return fIdx; }
    double cValue()   const { return fValue; }
    double cLow()     const { return lowEdge(fIdx); }
    double cHigh()    const { return highEdge(fIdx); }

    double getFullWidth()         const { return .5*(fMax + fMin); }
    double getFullCenter()        const { return (fMax - fMin); }
    double lowEdge  (int bin= 1)  const { return fMin+(bin-1)*(fMax-fMin)/fN; }
    double highEdge (int bin=-1)  const { if (bin==-1) bin=fN; return fMin+(bin)*(fMax-fMin)/fN; }
    double getCenter(int bin)     const { return (fMin + (bin-.5)*fW); }
    bool   isInside(double value) const { if(value>fMin && value<fMax) return true; return false; }

    void set     (double n, double min, double max, const char *ttl="");
    void setN    (double n)          { fN = n; fW = (fMax-fMin)/fN; }
    void setW    (double w)          { fW = w; fN = int((fMax-fMin)/fW); }
    void setMin  (double min)        { fMin = min; fW = (fMax-fMin)/fN; }
    void setMax  (double max)        { fMax = max; fW = (fMax-fMin)/fN; }
    void setTitle(const char *ttl)   { fTitle = ttl; }

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
      auto hist = new TH2D(nameHist,Form("%s;%s;%s;",titleHist,fTitleX,fTitleY),fNX,fMinX,fMaxX,fNY,fMinY,fMaxY);
      return hist;
    }
};



binning::binning(int n, double min, double max, double w) : fN(n), fMin(min), fMax(max), fW(w) { init(); }

binning::binning(TH1 *hist, int i)
{
  TAxis *axis;
       if (i==2) axis = hist -> GetYaxis();
  else if (i==3) axis = hist -> GetZaxis();
  else           axis = hist -> GetXaxis();

  fN = axis -> GetNbins();
  fMin = axis -> GetBinLowEdge(1);
  fMax = axis -> GetBinUpEdge(fN);
  fW = (fMax-fMin)/fN;
}

binning::binning(TGraph *graph, int i)
{
  fN = graph -> GetN();
  double x1,x2,y1,y2;
  graph -> ComputeRange(x1,y1,x2,y2);
  double xe1 = (x2-x1)/(fN-1)/2.;
  double xe2 = (x2-x1)/(fN-1)/2.;
  double ye1 = (y2-y1)/(fN-1)/2.;
  double ye2 = (y2-y1)/(fN-1)/2.;
  if (graph->InheritsFrom(TGraphErrors::Class())) {
    xe1 = graph -> GetErrorX(0);
    xe2 = graph -> GetErrorX(fN-1);
    ye1 = graph -> GetErrorY(0);
    ye2 = graph -> GetErrorY(fN-1);
  }
  if (i==2) { fMin = y1 - ye1; fMax = y2 + ye2; }
  else      { fMin = x1 - xe1; fMax = x2 + xe2; }
  fW = (fMax-fMin)/fN;
}

binning::binning(int n, double min, double max, const char *ttl) : fN(n), fMin(min), fMax(max), fTitle(ttl) { init(); }

void binning::reset(bool includeOUFlow)  { fIdx = (includeOUFlow?-1:0); }

bool binning::next(bool includeOUFlow) 
{
  if ((includeOUFlow && fIdx>fN) || (fIdx>fN-1))
    return false;
  fValue = fMin + (fIdx++) * fW + .5 * fW;
  return true;
}

TH1D *binning::newHist(const char *name, const char *title) {
  if (TString(title).IsNull()) title = fTitle;
  return (new TH1D(name,title,fN,fMin,fMax));
}

void binning::init()
{
       if (fW>0&&fN<1) setW(fW);
  else if (fN>0&&fW<1) setN(fN);
}


void binning::set(double n, double min, double max, const char *ttl)
{
  fN = n;
  fMin = min;
  fMax = max;
  fW = (fMax-fMin)/fN;
  if (TString(ttl).IsNull()==false)
    fTitle = ttl;
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
}

binning2 binning::operator*(const binning binn) {
  return binning2(fTitle,fN,fMin,fMax,binn.getTitle(),binn.getN(),binn.getMin(),binn.getMax());
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
