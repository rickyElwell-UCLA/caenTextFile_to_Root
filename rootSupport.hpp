

//////////////////////////////////////////////////////////
/// Model functions
//////////////////////////////////////////////////////////

Double_t fNoise(Double_t *x, Double_t *par);
Double_t fOnePeak(Double_t *x, Double_t *par);
Double_t fOnePeakClipped(Double_t *x, Double_t *par);
Double_t fTwoPeak(Double_t *x, Double_t *par);



//////////////////////////////////////////////////////////
/// Filtering Procedures
//////////////////////////////////////////////////////////

void findEventPairs(TTree *dev, TTree *pos_ch, TTree *neg_ch);
void iatByPosition(TTree *dev, TTree *iatByPosition);
