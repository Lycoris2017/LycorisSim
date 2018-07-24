
//double mycircle(std::vector<double>* x,std::vector<double>* y)
void Test_CircleFit(){
  
  std::vector<double>* x = new std::vector<double>;
  std::vector<double>* y = new std::vector<double>;

  x->push_back(2.);   y->push_back(0.); 
  x->push_back(0.);   y->push_back(2.); 
  x->push_back(1.);   y->push_back(sqrt(3.)); 
  //  x->push_back(0.0);   y->push_back(0.); 
  //  x->push_back(0.0);   y->push_back(5.); 
  //  x->push_back(1.0);   y->push_back(6.); 
  //  x->push_back(2.0);   y->push_back(7.); 
  //  x->push_back(8.0);   y->push_back(7.); 
  
  double * xArray = new double [x->size()];
  double * yArray = new double [x->size()];
  
  double * ui = new double [x->size()];
  double * vi = new double [x->size()];
  double * ui2 = new double [x->size()];
  double * vi2 = new double [x->size()];
  double * ui3 = new double [x->size()];
  double * vi3 = new double [x->size()];
  double * uv = new double [x->size()];
  double * uvv = new double [x->size()];
  double * vuu = new double [x->size()];

  double * Ri = new double [x->size()];
  
  for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
    {
      xArray[vectorIterator]=x->at(vectorIterator);
      yArray[vectorIterator]=y->at(vectorIterator);
    }
  
  double xMean=0., yMean=0.;
  xMean=TMath::Mean(x->size(),xArray);
  yMean=TMath::Mean(x->size(),yArray);
  
  //  cout<<" xmean "<<xMean<<" yMean "<<yMean<<endl;
  
  for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
    {
      ui[vectorIterator]=x->at(vectorIterator)-xMean;
      vi[vectorIterator]=y->at(vectorIterator)-yMean;
    }
  
  for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
    {
      ui2[vectorIterator]=ui[vectorIterator]*ui[vectorIterator];
      vi2[vectorIterator]=vi[vectorIterator]*vi[vectorIterator];
      ui3[vectorIterator]=ui[vectorIterator]*ui[vectorIterator]*ui[vectorIterator];
      vi3[vectorIterator]=vi[vectorIterator]*vi[vectorIterator]*vi[vectorIterator];
      
      uv[vectorIterator]=ui[vectorIterator]*vi[vectorIterator];
      uvv[vectorIterator]=ui[vectorIterator]*vi[vectorIterator]*vi[vectorIterator];
      vuu[vectorIterator]=vi[vectorIterator]*ui[vectorIterator]*ui[vectorIterator];
    }
  
  double Su=0., Suu=0., Suuu=0., Suv=0., Suvv=0.;
  double Sv=0., Svv=0., Svvv=0., Svuu=0.;
  
  for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
    {
      Su   += ui[vectorIterator];
      Sv   += vi[vectorIterator];
      Suu  += ui2[vectorIterator];
      Svv  += vi2[vectorIterator];
      Suuu += ui3[vectorIterator];
      Svvv += vi3[vectorIterator];
      Suv  += uv[vectorIterator];
      Suvv += uvv[vectorIterator];
      Svuu += vuu[vectorIterator];
    }
  
  //cout<<"Suu="<<Suu<<", Suv="<<Suv<<", Svv="<<Svv<<", Suuu="<<Suuu<<", Svvv="<<Svvv<<", Suvv="<<Suvv<<", Svuu="<<Svuu<<endl; 
  
  double vc=0., uc=0., a=0.;
  
  vc = ( Suv*(Suuu+Suvv) - Suu*(Svvv+Svuu) ) / (2* (Suv*Suv - Svv*Suu) );
  uc = ( 0.5* (Suuu+Suvv) - vc*Suv) / Suu;
  
  //cout<<"uc "<<uc<<" vc "<<vc<<endl;

  // here we get the center of the circle (xc,yc)  
  double xc=0., yc=0.;
  xc= uc + xMean;
  yc = vc + yMean;
  
  int Nmeas = x->size();
  // and using also the number of measurements, we get back the radius sqrt(a)
  double N = double(Nmeas);
  a = vc*vc + uc*uc + (Suu+Svv/N); // need to define N
  
  cout<<" xcenter = "<<xc<<", y centre = "<<yc<<" radius = "<<sqrt(a)<<endl;


  for (int i=0; i<x->size(); i++ ) {
    Ri[i] = sqrt((x->at(i)-xc)**2 + (y->at(i)-yc)**2);
  }
  cout<<"ok1"<<endl;
  double Rcircle = TMath::Mean(x->size(),Ri);
  cout<<"ok2"<<endl;
  //  residu_1 = TMath::Sum((Ri-Rcircle)**2);
  cout<<"Ri "<<Ri[0]<<" "<<Ri[1]<<" "<<Ri[2]<<". Rcircle "<<Rcircle<<endl;

  // Now I have the center of the circle but I want to get the distance from the true hits in the TPC
  //first get vector with hits in the TPC
  std::vector<double>* xTPC = new std::vector<double>;
  std::vector<double>* yTPC = new std::vector<double>;

  xTPC->push_back(0.0);   yTPC->push_back(0.); 
  
  double * distance = new double [xTPC->size()];
  double * resid = new double [xTPC->size()];

  // then calculate distance from circle center
  //  double distance=0., resid=0.;
  double radius = sqrt(a);
  for (unsigned int vectorIterator=0;vectorIterator<xTPC->size();vectorIterator++) 
    {
      distance[vectorIterator] = sqrt( (xTPC->at(vectorIterator)-xc)**2 + (yTPC->at(vectorIterator)-yc)**2 );
      if ( distance[vectorIterator] > radius ) resid[vectorIterator] = distance[vectorIterator] - radius;
      else resid[vectorIterator] = radius - distance[vectorIterator];
    }

  cout<<distance[0]<<" "<<radius<<" "<<resid[0]<<endl;
}
