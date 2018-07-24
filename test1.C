double calculateMomentum(std::vector<double>* x,std::vector<double>* y)
{
	double * xArray = new double [x->size()];
	double * yArray = new double [x->size()];
	double * x2 = new double [x->size()];
	double * y2 = new double [x->size()];
	double * xy = new double [x->size()];
	double * xr2 = new double [x->size()];
	double * yr2 = new double [x->size()];
	double * r2 = new double [x->size()];
	double * r4 = new double [x->size()];
	for(unsigned int vectorIterator=0;vectorIterator<x->size();vectorIterator++)
	{
		xArray[vectorIterator]=x->at(vectorIterator);
		yArray[vectorIterator]=y->at(vectorIterator);
		x2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator);
		y2[vectorIterator]=y->at(vectorIterator)*y->at(vectorIterator);
		xy[vectorIterator]=x->at(vectorIterator)*y->at(vectorIterator);
		r2[vectorIterator]=x->at(vectorIterator)*x->at(vectorIterator)+y->at(vectorIterator)*y->at(vectorIterator);
		xr2[vectorIterator]=x->at(vectorIterator)*r2[vectorIterator];
		yr2[vectorIterator]=y->at(vectorIterator)*r2[vectorIterator];
		r4[vectorIterator]=r2[vectorIterator]*r2[vectorIterator];
	}
	double Cxx,Cxy,Cyy,Cxr2,Cyr2,Cr2r2,phi,kappa,delta,rho,momentum;
	double xMean,yMean,x2Mean,y2Mean,xyMean,xr2Mean,yr2Mean,r2Mean,r4Mean;
	xMean=TMath::Mean(x->size(),xArray);

	x2Mean=TMath::Mean(x->size(),x2);
	yMean=TMath::Mean(x->size(),yArray);
	y2Mean=TMath::Mean(x->size(),y2);
	xyMean=TMath::Mean(x->size(),xy);
	xr2Mean=TMath::Mean(x->size(),xr2);
	yr2Mean=TMath::Mean(x->size(),yr2);
	r2Mean=TMath::Mean(x->size(),r2);
	r4Mean=TMath::Mean(x->size(),r4);
	Cxx=x2Mean-(xMean*xMean);
	Cxy=xyMean-xMean*yMean;
	Cyy=y2Mean-yMean*yMean;
	Cxr2=xr2Mean-xMean*r2Mean;
	Cyr2=yr2Mean-yMean*r2Mean;
	Cr2r2=r4Mean-r2Mean*r2Mean;

	phi=0.5*TMath::ATan(2*(Cr2r2*Cxy-Cxr2*Cyr2)/(Cr2r2*(Cxx-Cyy)-Cxr2*Cxr2+Cyr2*Cyr2));
	kappa=(TMath::Sin(phi)*Cxr2-TMath::Cos(phi)*Cyr2)/Cr2r2;
	delta=-kappa*r2Mean+TMath::Sin(phi)*xMean-TMath::Cos(phi)*yMean;

	rho=2*kappa/(TMath::Sqrt(1-4*delta*kappa));
	momentum=0.3*1/rho;
	delete [] xArray;
	delete [] yArray;
	delete [] x2;
	delete [] y2;
	delete [] xy;
	delete [] xr2;
	delete [] yr2;
	delete [] r2;
	delete [] r4;
	return (momentum);
};
