//double calculateMomentum(std::vector<double>* x,std::vector<double>* y)
//{
//	double * xArray = new double [x->size()];
//	double * yArray = new double [x->size()];
//	double * x2 = new double [x->size()];
//	double * y2 = new double [x->size()];
//	double * xy = new double [x->size()];
//	double * xr2 = new double [x->size()];
//	double * yr2 = new double [x->size()];
//	double * r2 = new double [x->size()];
//	double * r4 = new double [x->size()];
//	for(int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
//	{
//		//std::cout<<"test0"<<std::endl;
//		xArray[vectorIterator]=x->at(vectorIterator);
//		//std::cout<<"test1"<<std::endl;
//		yArray[vectorIterator]=y->at(vectorIterator);
//		x2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator);
//		y2[vectorIterator]=y->at(vectorIterator)*y->at(vectorIterator);
//		xy[vectorIterator]=x->at(vectorIterator)*y->at(vectorIterator);
//		r2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator)+y->at(vectorIterator)*y->at(vectorIterator);
//		xr2[vectorIterator]=x->at(vectorIterator)*r2[vectorIterator];
//		yr2[vectorIterator]=y->at(vectorIterator)*r2[vectorIterator];
//		r4[vectorIterator]=r2[vectorIterator]*r2[vectorIterator];
//	}
//	double Cxx,Cxy,Cyy,Cxr2,Cyr2,Cr2r2,phi,kappa,delta,rho,momentum;
//	double xMean,yMean,x2Mean,xMean,y2Mean,xyMean,xr2Mean,yr2Mean,r2Mean,r4Mean;
//	xMean=TMath::Mean(x->size(),xArray);
//
//	x2Mean=TMath::Mean(x->size(),x2);
//	yMean=TMath::Mean(x->size(),yArray);
//	y2Mean=TMath::Mean(x->size(),y2);
//	xyMean=TMath::Mean(x->size(),xy);
//	xr2Mean=TMath::Mean(x->size(),xr2);
//	yr2Mean=TMath::Mean(x->size(),yr2);
//	r2Mean=TMath::Mean(x->size(),r2);
//	r4Mean=TMath::Mean(x->size(),r4);
//	Cxx=x2Mean-(xMean*xMean);
//	Cxy=xyMean-xMean*yMean;
//	Cyy=y2Mean-yMean*yMean;
//	Cxr2=xr2Mean-xMean*r2Mean;
//	Cyr2=yr2Mean-yMean*r2Mean;
//	Cr2r2=r4Mean-r2Mean*r2Mean;
//
//	phi=0.5*TMath::ATan(2*(Cr2r2*Cxy-Cxr2*Cyr2)/(Cr2r2*(Cxx-Cyy)-Cxr2*Cxr2+Cyr2*Cyr2));
//	kappa=(TMath::Sin(phi)*Cxr2-TMath::Cos(phi)*Cyr2)/Cr2r2;
//	delta=-kappa*r2Mean+TMath::Sin(phi)*xMean-TMath::Cos(phi)*yMean;
//
//	rho=2*kappa/(TMath::Sqrt(1-4*delta*kappa));
//	momentum=0.3*1/rho;
//	delete [] xArray;
//	delete [] yArray;
//	delete [] x2;
//	delete [] y2;
//	delete [] xy;
//	delete [] xr2;
//	delete [] yr2;
//	delete [] r2;
//	delete [] r4;
//	return (momentum);
//};
//
////function to calculate the momentum of a track with weighted points. You just need to insert the x and y positions and x and y error as std::vectors in units of mm to get the momentum in MeV
//double calculateWeightedMomentum(std::vector<double>* x,std::vector<double>* y,std::vector<double>* xerr,std::vector<double>* yerr)
//{
//	double * xArray = new double [x->size()];
//	double * yArray = new double [x->size()];
//	double * x2 = new double [x->size()];
//	double * y2 = new double [x->size()];
//	double * xy = new double [x->size()];
//	double * r2 = new double [x->size()];
//	double * xr2 = new double [x->size()];
//	double * yr2 = new double [x->size()];
//	double * r4 = new double [x->size()];
//	double * weight = new double [x->size()];
//	for(int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
//	{
//		//std::cout<<"test0"<<std::endl;
//		weight[vectorIterator]=1/(yerr->at(vectorIterator)*yerr->at(vectorIterator));
//		xArray[vectorIterator]=x->at(vectorIterator);
//		yArray[vectorIterator]=y->at(vectorIterator);
//		x2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator);
//		y2[vectorIterator]=y->at(vectorIterator)*y->at(vectorIterator);
//		xy[vectorIterator]=x->at(vectorIterator)*y->at(vectorIterator);
//		r2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator)+y->at(vectorIterator)*y->at(vectorIterator);
//		xr2[vectorIterator]=x->at(vectorIterator)*r2[vectorIterator];
//		yr2[vectorIterator]=y->at(vectorIterator)*r2[vectorIterator];
//		r4[vectorIterator]=r2[vectorIterator]*r2[vectorIterator];
//	}
//	double Cxx,Cxy,Cyy,Cxr2,Cyr2,Cr2r2,phi,kappa,delta,rho,momentum;
//	double xMean,yMean,x2Mean,xMean,y2Mean,xyMean,xr2Mean,yr2Mean,r2Mean,r4Mean;
//	xMean=TMath::Mean(x->size(),xArray,weight);
//	x2Mean=TMath::Mean(x->size(),x2,weight);
//	yMean=TMath::Mean(x->size(),yArray,weight);
//	y2Mean=TMath::Mean(x->size(),y2,weight);
//	xyMean=TMath::Mean(x->size(),xy,weight);
//	xr2Mean=TMath::Mean(x->size(),xr2,weight);
//	yr2Mean=TMath::Mean(x->size(),yr2,weight);
//	r2Mean=TMath::Mean(x->size(),r2,weight);
//	r4Mean=TMath::Mean(x->size(),r4,weight);
//
//	Cxx=x2Mean-(xMean*xMean);
//	Cxy=xyMean-xMean*yMean;
//	Cyy=y2Mean-yMean*yMean;
//	Cxr2=xr2Mean-xMean*r2Mean;
//	Cyr2=yr2Mean-yMean*r2Mean;
//	Cr2r2=r4Mean-r2Mean*r2Mean;
//
//	phi=0.5*TMath::ATan(2*(Cr2r2*Cxy-Cxr2*Cyr2)/(Cr2r2*(Cxx-Cyy)-Cxr2*Cxr2+Cyr2*Cyr2));
//	kappa=(TMath::Sin(phi)*Cxr2-TMath::Cos(phi)*Cyr2)/Cr2r2;
//	delta=-kappa*r2Mean+TMath::Sin(phi)*xMean-TMath::Cos(phi)*yMean;
//
//	rho=2*kappa/(TMath::Sqrt(1-4*delta*kappa));
//	momentum=0.3*1/rho;
//	delete [] xArray;
//	delete [] weight;
//	delete [] yArray;
//	delete [] x2;
//	delete [] y2;
//	delete [] xy;
//	delete [] xr2;
//	delete [] yr2;
//	delete [] r2;
//	delete [] r4;
//
//	return (momentum);
//};
