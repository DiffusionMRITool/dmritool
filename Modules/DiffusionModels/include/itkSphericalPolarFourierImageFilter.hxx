/**
 *       @file  itkSphericalPolarFourierImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-26-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkSphericalPolarFourierImageFilter_hxx
#define __itkSphericalPolarFourierImageFilter_hxx


#include "itkSphericalPolarFourierImageFilter.h"
#include "itkSPFScaleFromMeanDiffusivityImageFilter.h"
#include "utlVNLBlas.h"
#include "utl.h"

#include <itkProgressReporter.h>

namespace itk
{

template< class TInputImage, class TOutputImage >
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::SphericalPolarFourierImageFilter() : Superclass(), 
  m_Gn0(new STDVectorType()),
  m_G0DWI(new VectorType())
{
  // NOTE: ComputeScale is a virtual function, thus it must be called in derived calss, not in base class
  this->ComputeScale(true);
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::SetBasisScale(const double scale)
{
  itkShowPositionThreadedLogger(this->GetDebug());
  double scale_old = this->m_BasisScale;  
  if (scale>0)
    this->m_BasisScale = scale;
  else
    this->ComputeScale(true);
  itkDebugMacro("setting this->m_BasisScale to " << this->m_BasisScale);

  if (scale>0 && std::fabs((scale_old-this->m_BasisScale)/this->m_BasisScale)>1e-8)
    {
    this->Modified();
    this->m_BasisRadialMatrix=MatrixPointer(new MatrixType());
    this->m_BasisMatrix=MatrixPointer(new MatrixType());
    this->m_BasisMatrixForB0=MatrixPointer(new MatrixType());
    this->m_Gn0=STDVectorPointer(new STDVectorType());
    this->m_G0DWI=VectorPointer(new VectorType());
    }
}

template< class TInputImage, class TOutputImage >
typename LightObject::Pointer
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::InternalClone() const
{
  itkShowPositionThreadedLogger(this->GetDebug());
  typename LightObject::Pointer loPtr = Superclass::InternalClone();

  typename Self::Pointer rval = dynamic_cast<Self *>(loPtr.GetPointer());
  if(rval.IsNull())
    {
    itkExceptionMacro(<< "downcast to type " << this->GetNameOfClass()<< " failed.");
    }
  
  rval->m_Gn0 = m_Gn0;
  rval->m_G0DWI = m_G0DWI;
  return loPtr;
}

template< class TInputImage, class TOutputImage >
double
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeScale(const bool setScale)
{
  itkShowPositionThreadedLogger(this->GetDebug());
  double scale;
  double tau = this->m_SamplingSchemeQSpace->GetTau();
  if (this->m_IsOriginalBasis)
    scale = 1.0 / (8*M_PI*M_PI*tau*this->m_MD0);  // 714.29 (700) scale for SPF basis, dual scale is 1/(4*pi^2*scale)
  else
    scale = 2*tau*this->m_MD0;  // 3.5462e-5 scale for SPF basis, dual scale is 1/(4*pi^2*scale)

  // std::cout << "scale=" << scale << std::endl << std::flush;
  // utlPrintVar3(true, this->m_IsOriginalBasis, tau, this->m_MD0);
  if (setScale)
    this->SetBasisScale(scale);
  if (this->GetDebug())
    std::cout << "m_BasisScale = " << this->m_BasisScale << std::endl;
  return scale;
}

template< class TInputImage, class TOutputImage >
std::vector<int>
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::GetIndexNLM ( const int index ) const
{
  int sh_b = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n = index / sh_b;
  int residual = index - n*sh_b;
  int j=0;
  std::vector<int> nlm;
  nlm.push_back(n);
  std::vector<int> lm = utl::GetIndexSHlm(residual);
  nlm.push_back(lm[0]);
  nlm.push_back(lm[1]);
  return nlm;
}

template< class TInputImage, class TOutputImage >
int
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::GetIndexJ (const int n, const int l, const int m) const
{
  return n*(this->m_SHRank+1)*(this->m_SHRank+2)/2 + utl::GetIndexSHj(l,m);
}

template< class TInputImage, class TOutputImage >
std::vector<int>
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::DimToRank ( const int dimm ) const
{
  std::vector<int> result;
  int radialRank=-1, shRank=-1;
  for ( int radialRank = 0; radialRank <= 10; radialRank += 1 ) 
    {
    for ( int shRank = 0; shRank <= 12; shRank += 2 ) 
      {
      int dim = RankToDim(false, radialRank, shRank);
      if (dim==dimm)
        {
        result.push_back(radialRank);
        result.push_back(shRank);
        return result;
        }
      }
    }
  utlException(true, "wrong logic");
  return result;
}   

template< class TInputImage, class TOutputImage >
int
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::RankToDim (const bool is_radial, const int radialRank, const int shRank) const
{
  int radialRank_real = radialRank>=0?radialRank:this->m_RadialRank;
  int shRank_real = shRank>=0?shRank:this->m_SHRank;
  utlException(radialRank_real<0 || shRank_real<0, "wrong rank");
  if (is_radial)
    return radialRank_real+1;
  else
    return (shRank_real + 1)*(shRank_real + 2)/2*(radialRank_real+1);
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeRadialVectorForE0InBasis ()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  typename SPFGenerator::Pointer spf = SPFGenerator::New();
  m_Gn0 = STDVectorPointer(new STDVectorType(this->m_RadialRank+1,-100));
  spf->SetSPFType(SPFGenerator::SPF);
  spf->SetScale(this->m_BasisScale);
  for ( int i = 0; i < this->m_RadialRank+1; i += 1 ) 
    {
    spf->SetN(i);
    (*m_Gn0)[i] = spf->Evaluate(0);
    }
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeRadialVectorForE0InDWI ()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  if (this->m_SamplingSchemeQSpace->GetIndicesInShells()->size()==0)
    this->m_SamplingSchemeQSpace->GroupRadiusValues();

  m_G0DWI = VectorPointer(new VectorType(this->m_SamplingSchemeQSpace->GetNumberOfSamples(),-100));
  typename Superclass::SamplingSchemeQSpaceType::Index2DVectorPointer indices = this->m_SamplingSchemeQSpace->GetIndicesInShells();
  STDVectorPointer qVector = this->m_SamplingSchemeQSpace->GetRadiusVector();
  for ( int i = 0; i < indices->size(); i += 1 ) 
    {
    typename Superclass::SamplingSchemeQSpaceType::IndexVectorType indexTemp = (*indices)[i];
    double qq = (*qVector)[ indexTemp[0] ];
    double val = std::exp(-qq*qq/(2.0*this->m_BasisScale));
    for ( int j = 0; j < indexTemp.size(); j += 1 ) 
      (*m_G0DWI)[ indexTemp[j] ] = val;
    }
}

// template< class TInputImage, class TOutputImage >
// void
// SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
// ::VerifyInputParameters() const
// {
//   itkShowPositionThreadedLogger(this->GetDebug());
//   Superclass::VerifyInputParameters();
// }

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeRadialMatrix ()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  
  if (this->m_SamplingSchemeQSpace->GetBVector()->size()>0 && this->m_SamplingSchemeQSpace->GetRadiusVector()->size()==0)
    this->m_SamplingSchemeQSpace->ConvertBVectorToQVector();
  
  if (this->m_SamplingSchemeQSpace->GetIndicesInShells()->size()==0)
    this->m_SamplingSchemeQSpace->GroupRadiusValues();
  
  // NOTE: m_BVector.size()==m_QOrientations.size()
  int n_s, n_b;
  n_s = this->m_SamplingSchemeQSpace->GetRadiusVector()->size();
  utlGlobalException( n_s==0, "no q vector");
  
  const STDVectorPointer qVector = this->m_SamplingSchemeQSpace->GetRadiusVector();
  typename Superclass::SamplingSchemeQSpaceType::Index2DVectorPointer indices = this->m_SamplingSchemeQSpace->GetIndicesInShells();

  std::string threadIDStr = this->ThreadIDToString();
  if (this->GetDebug())
    {
    std::ostringstream msg;
    msg << threadIDStr << "m_BasisScale = " << this->m_BasisScale << ", m_Tau = " << this->m_SamplingSchemeQSpace->GetTau() << ", numberOfShell = " << indices->size() << std::endl;
    this->WriteLogger(msg.str());
    }
    
  MatrixPointer B;
  typename SPFGenerator::Pointer spf = SPFGenerator::New();
  spf->SetScale(this->m_BasisScale);

  if (this->m_IsOriginalBasis)
    {
    spf->SetSPFType(SPFGenerator::SPF);

    // the basis of order N has N+1 terms (0,1,...,N)
    n_b = this->m_RadialRank + 1;


    if(this->GetDebug())
      {
      std::ostringstream msg;
      msg << threadIDStr << "Generating the "<< n_s << "x" << n_b << " RadialMatrix...\n";
      this->WriteLogger(msg.str());
      }

    B = MatrixPointer(new MatrixType(n_s,n_b));

    for ( int ib = 0; ib < n_b; ib += 1 ) 
      {
      spf->SetN(ib);
      for ( int shell = 0; shell < indices->size(); shell += 1 ) 
        {
        typename Superclass::SamplingSchemeQSpaceType::IndexVectorType indexTemp = (*indices)[shell];
        double qq = (*qVector)[ indexTemp[0] ];
        double spfVal = spf->Evaluate(qq,false);
        for ( int js = 0; js < indexTemp.size(); js += 1 ) 
          (*B)(indexTemp[js],ib) = spfVal;
        }
      }

    }
  else
    {
    spf->SetSPFType(SPFGenerator::DSPF);

    n_b = (this->m_RadialRank+1)*(this->m_SHRank/2+1);

    if(this->GetDebug())
      {
      std::ostringstream msg;
      msg << threadIDStr << "Generating the "<< n_s << "x" << n_b << " RadialMatrix...\n";
      this->WriteLogger(msg.str());
      }

    B = MatrixPointer(new MatrixType(n_s,n_b));

    for ( int ib = 0; ib < this->m_RadialRank+1; ib += 1 ) 
      {
      spf->SetN(ib);
      for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
        {
        spf->SetL(l);
        int col = ib*(this->m_SHRank/2+1)+l/2;
        for ( int shell = 0; shell < indices->size(); shell += 1 ) 
          {
          typename Superclass::SamplingSchemeQSpaceType::IndexVectorType indexTemp = (*indices)[shell];
          double qq = (*qVector)[ indexTemp[0] ];
          double spfVal = spf->Evaluate(qq,false);
          for ( int js = 0; js < indexTemp.size(); js += 1 ) 
            (*B)(indexTemp[js],col) = spfVal;
          }
        }
      }

    }
  
  if(this->GetDebug()) 
    {
    std::ostringstream msg;
    utl::PrintUtlMatrix(*B,"RadialMatrix", " ", msg << threadIDStr);
    this->WriteLogger(msg.str());
    }

  this->m_BasisRadialMatrix = B;
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeBasisMatrix() 
{
  itkShowPositionThreadedLogger(this->GetDebug());

  if (this->m_SamplingSchemeQSpace->GetBVector()->size()>0 && this->m_SamplingSchemeQSpace->GetRadiusVector()->size()==0)
    this->m_SamplingSchemeQSpace->ConvertBVectorToQVector();
  
  if (this->m_SamplingSchemeQSpace->GetIndicesInShells()->size()==0)
    this->m_SamplingSchemeQSpace->GroupRadiusValues();

  if (this->m_BasisSHMatrix->Rows()==0)
    this->ComputeSHMatrix();
  MatrixPointer B_sh = this->m_BasisSHMatrix;

  if (this->m_BasisRadialMatrix->Rows()==0)
    this->ComputeRadialMatrix();
  MatrixPointer B_ra = this->m_BasisRadialMatrix;

  MatrixPointer qOrientations = this->m_SamplingSchemeQSpace->GetOrientationsSpherical();

  int qVector_size = this->m_SamplingSchemeQSpace->GetRadiusVector()->size();
  int grad_size = qOrientations->Rows();
  
  std::string threadIDStr = this->ThreadIDToString();
  if (this->GetDebug())
    {
    std::ostringstream msg;
    msg << threadIDStr << "this->m_SamplingSchemeQSpace->GetRadiusVector()->size() = " << this->m_SamplingSchemeQSpace->GetRadiusVector()->size() << std::endl;
    msg << threadIDStr << "qOrientations->Rows() = " << qOrientations->Rows() << std::endl;
    msg << threadIDStr << "qVector_size = " << qVector_size << std::endl;
    msg << threadIDStr << "grad_size = " << grad_size << std::endl;
    this->WriteLogger(msg.str());
    }
  utlException(qVector_size!=grad_size, "qVector_size and grad_size should keep the same size");
  int n_s = qVector_size;

  // the basis has two parts 
  int n_b_ra = this->m_RadialRank + 1;
  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n_b =  n_b_ra * n_b_sh;

  if(this->GetDebug())
    {
    std::ostringstream msg;
    msg << threadIDStr << "n_b_ra = " << n_b_ra  << std::endl;
    msg << threadIDStr << "n_b_sh = " << n_b_sh  << std::endl;
    msg << threadIDStr << "n_b = " << n_b  << std::endl;
    msg << threadIDStr << "Generating the "<< n_s << "x" << n_b << " basis matrix...\n";
    this->WriteLogger(msg.str());
    }

  MatrixPointer B(new MatrixType(n_s, n_b));
  utlException(B_sh->Columns()!=n_b_sh, "the SHMatrix does not have the right width. B_sh=" << *B_sh << ", n_b_sh="<< n_b_sh);
  utlException(B_sh->Rows()!=B_ra->Rows(), "the SHMatrix and the RadialMatrix do not have the same samples. B_sh="<< *B_sh << ", B_ra=" << *B_ra);
  // B_sh.print("B_sh",2);
  // B.print("B",2);
  utlException(B_sh->Rows()!=B->Rows(), "the SHMatrix and the basisMatrix do not have the same samples. B_sh="<< *B_sh << ", B="<< *B);

  if (this->m_IsOriginalBasis)
    {
    utlException(B_ra->Columns()!=n_b_ra, "the RadialMatrix does not have the right width");

    // utl::Tic(std::cout<<"index start");
    double *B_sh_data = B_sh->GetData();
    double *B_ra_data = B_ra->GetData();
    double *B_data = B->GetData();
    int index_ra=0, index_B=0;
    std::vector<int> index_sh(n_b_ra,0);
    for ( int k = 0; k < n_s; k += 1 ) 
      {
      for ( int i = 0; i < n_b_ra; i += 1 ) 
        {
        for ( int j = 0; j < n_b_sh; j += 1 ) 
          {
          B_data[index_B] = B_ra_data[index_ra] * B_sh_data[index_sh[i] ]; 
          // (*B)(k,i*n_b_sh+j) = (*B_ra_data)(k,i) * (*B_sh)(k,j); 
          index_sh[i]++;
          index_B++;
          }
        index_ra++;
        }
      }
    // utl::Toc();

    // utl::Tic(std::cout<<"() start");
    // for ( int i = 0; i < n_b_ra; i += 1 ) 
    //   for ( int j = 0; j < n_b_sh; j += 1 ) 
    //     for ( int k = 0; k < n_s; k += 1 ) 
    //       {
    //       (*B)(k,i*n_b_sh+j) = (*B_ra)(k,i) * (*B_sh)(k,j); 
    //       }
    // utl::Toc();
    }
  else
    {
    utlException(B_ra->Columns()!=n_b_ra*(this->m_SHRank/2+1), "the RadialMatrix does not have the right width");
    for ( int n = 0; n < n_b_ra; n += 1 ) 
      {
      int jj = 0;
      for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
        {
        int col_ra = n*(this->m_SHRank/2+1)+l/2;
        for ( int m = -l; m <= l; m += 1 ) 
          {
          int col = n*n_b_sh+jj;
          for ( int k = 0; k < n_s; k += 1 ) 
            {
            (*B)(k,col) = (*B_ra)(k,col_ra) * (*B_sh)(k,jj); 
            }
          jj++;
          }
        }
      }
    }

  // // add large weight for constrains b=0,  (now bVector has no 0)
  // if ( !this->m_IsAnalyticalB0 && std::abs(this->m_B0Weight-1)>1e-8 && this->m_B0Weight>0)
  //   {
  //   std::vector<int> index = utl::FindVector(*this->m_SamplingSchemeQSpace->GetRadiusVector() , double(0) , double(1e-6));
  //   if (this->GetDebug())
  //     {
  //     std::cout << "use this->m_B0Weight = " << this->m_B0Weight << "when computing basis matrix" << std::endl;
  //     utl::PrintVector(index,"index");
  //     }
  //   for ( int i = 0; i < index.size(); i += 1 ) 
  //     {
  //     for ( int j = 0; j < B->Columns(); j += 1 ) 
  //       {
  //       (*B)(index[i],j) *= this->m_B0Weight;
  //       }
  //     }
  //   }
  
  if(this->GetDebug()) 
    {
    std::ostringstream msg;
    utl::PrintUtlMatrix(*B,"BasisMatrix", " ", msg<< threadIDStr);
    this->WriteLogger(msg.str());
    }

  this->m_BasisMatrix = B;
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeBasisMatrixForB0 ()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  Pointer selfClone = this->Clone();
  selfClone->m_BasisMatrix = MatrixPointer(new MatrixType());
  selfClone->m_BasisRadialMatrix = MatrixPointer(new MatrixType());
  selfClone->m_BasisSHMatrix = MatrixPointer(new MatrixType());
  selfClone->m_SamplingSchemeQSpace = SamplingSchemeQSpaceType::New();
  SamplingSchemeQSpacePointer sampling = selfClone->GetSamplingSchemeQSpace();
  static MatrixPointer grad = utl::ReadGrad<double>(3, DIRECTION_NODUPLICATE, CARTESIAN_TO_CARTESIAN);
  STDVectorPointer b0Vector = STDVectorPointer( new STDVectorType(grad->Rows(),0.0));
  sampling->SetBVector(b0Vector);
  sampling->SetOrientationsCartesian(grad);
  selfClone->ComputeBasisMatrix();
  this->m_BasisMatrixForB0 = selfClone->GetBasisMatrix();
}    // -----  end of method EAP<T>::<method>  -----

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ComputeRegularizationWeight( )
{
  itkShowPositionThreadedLogger(this->GetDebug());

  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n_b_ra = this->m_RadialRank + 1;

  utlException((this->m_EstimationType==Self::LS || this->m_EstimationType==Self::L1_2) && this->m_LambdaSpherical<0 && this->m_LambdaRadial<0, "need to set m_LambdaSpherical and m_LambdaSpherical");
  utlException((this->m_EstimationType==Self::L1_DL) && this->m_LambdaL1<0, "need to set m_LambdaL1");
  if (this->m_EstimationType!=Self::L1_DL)
    {
    this->m_RegularizationWeight=VectorPointer( new VectorType(!this->m_IsAnalyticalB0?this->RankToDim():(n_b_sh*this->m_RadialRank)) );
    for ( int i = 0; i <= ((!this->m_IsAnalyticalB0)?this->m_RadialRank:(this->m_RadialRank-1)); i += 1 ) 
      {
      int j = 0;
      for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
        {
        double lambda=0;
        if (this->m_IsAnalyticalB0)
          lambda = this->m_LambdaSpherical*l*l*(l+1)*(l+1) + this->m_LambdaRadial*(i+1)*(i+1)*(i+2)*(i+2);  // use order 1,2,...,N, do not use i=0
        else
          lambda = this->m_LambdaSpherical*l*l*(l+1)*(l+1) + this->m_LambdaRadial*i*i*(i+1)*(i+1);
        for ( int m = -l; m <= l; m += 1 ) 
          {
          if (this->m_EstimationType==Self::L1_2 || this->m_EstimationType==Self::LS)
            (*this->m_RegularizationWeight)[i*n_b_sh+j] = lambda;
          j++;
          }
        }
      }
    }
  else
    {
    if (this->m_BasisCombinationMatrix->Size()==0)
      {
      utl::ReadMatrix(utl::CreateExpandedPath(utl::LearnedSPFDictionary_SH8_RA4_K250), *this->m_BasisCombinationMatrix);
      }
    if (this->GetDebug())
      utl::PrintUtlMatrix(*this->m_BasisCombinationMatrix, "*this->m_BasisCombinationMatrix");

    if (this->m_BasisEnergyDL->Size()==0)
      {
      if (std::abs(this->m_BasisEnergyPowerDL)>1e-10)
        {
        std::vector<double> vecTemp;
        utl::ReadVector(utl::CreateExpandedPath(utl::LearnedSPFEnergy_SH8_RA4_K250), vecTemp);
        *this->m_BasisEnergyDL = utl::StdVectorToUtlVector(vecTemp);
        *this->m_BasisEnergyDL /= this->m_BasisEnergyDL->GetMean();

        if (std::abs(this->m_BasisEnergyPowerDL-1)>1e-10)
          utl::PowerVector(this->m_BasisEnergyDL->Begin(), this->m_BasisEnergyDL->End(), this->m_BasisEnergyPowerDL);
        }
      else
        {
        this->m_BasisEnergyDL=VectorPointer(new VectorType(this->m_BasisCombinationMatrix->Columns()));
        this->m_BasisEnergyDL->Fill(1.0);
        }
      }
    if (this->GetDebug())
      utl::PrintUtlVector(*this->m_BasisEnergyDL, "*this->m_BasisEnergyDL");

    utlException(this->m_BasisMatrix->Size()>0 && this->m_BasisCombinationMatrix->Rows()!=this->m_BasisMatrix->Columns()-n_b_sh, "wrong size of dictionary. m_BasisCombinationMatrix->Rows()="<<this->m_BasisCombinationMatrix->Rows() <<", BasisMatrix->Columns()-n_b_sh="<< this->m_BasisMatrix->Columns()-n_b_sh);
    utlException(this->m_BasisMatrix->Size()>0 && this->m_BasisCombinationMatrix->Columns()!=this->m_BasisEnergyDL->Size(), "wrong size of dictionary. this->m_BasisEnergyDL.size()="<<this->m_BasisEnergyDL->Size());

    this->m_RegularizationWeight=VectorPointer(new VectorType(this->m_BasisEnergyDL->Size()) );
    for ( int i = 0; i < this->m_BasisEnergyDL->Size(); i += 1 ) 
      {
      if (this->m_EstimationType==Self::L1_DL)
        (*this->m_RegularizationWeight)[i] = this->m_LambdaL1 / (*this->m_BasisEnergyDL)[i];
      else 
        (*this->m_RegularizationWeight)[i] = this->m_LambdaL2 / (*this->m_BasisEnergyDL)[i];
      }
    }

}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  itkShowPositionThreadedLogger(this->GetDebug());
  Superclass::BeforeThreadedGenerateData();

  if (this->m_SamplingSchemeQSpace->GetBVector()->size()>0 && this->m_SamplingSchemeQSpace->GetRadiusVector()->size()==0)
    this->m_SamplingSchemeQSpace->ConvertBVectorToQVector();

  if (this->m_SamplingSchemeQSpace->GetIndicesInShells()->size()==0)
    this->m_SamplingSchemeQSpace->GroupRadiusValues();

  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n_b_ra = this->m_RadialRank + 1;
  typename SPFGenerator::Pointer spf = SPFGenerator::New();

  if (!this->m_IsAnalyticalB0)
    {
    std::cout << "Use numerical way for constraint E(0)=1" << std::endl << std::flush;
    this->ComputeBasisMatrixForB0();
    }
  else
    {
    if (m_Gn0->size()==0)
      {
      this->ComputeRadialVectorForE0InBasis();
      this->ComputeRadialVectorForE0InDWI();
      }
    utlGlobalException(!this->m_IsOriginalBasis, "TODO: use analytical way for E(0)=1 and DSPF basis");
    std::cout << "Use analytical way for constraint E(0)=1" << std::endl << std::flush;
    }

  if (!this->IsAdaptiveScale())
    {
    std::cout << "Use the same scale for all voxels!" << std::endl << std::flush;
    this->ComputeBasisMatrix();

    MatrixPointer basisMatrix(new MatrixType());
    if (!this->m_IsAnalyticalB0)
      {
      utlGlobalException(this->m_EstimationType==Self::L1_DL, "L1-DL only supports analytical way");
      *basisMatrix = utl::ConnectUtlMatrix(*this->m_BasisMatrix, *utl::ToMatrix<double>(*this->m_BasisMatrixForB0 % this->m_B0Weight), true);
      }
    else
      {
      basisMatrix = MatrixPointer(new MatrixType(this->m_BasisMatrix->Rows(), n_b_sh*(n_b_ra-1)));
      utlException((*m_Gn0)[0]==0, "it should be not zero!, (*m_Gn0)[0]="<< (*m_Gn0)[0]);
      // utlPrintVar2 (this->_is_b0_analytical, this->m_RadialRank, R_N_0);

      for ( int i = 0; i < this->m_RadialRank; i += 1 ) 
        {
        for ( int j = 0; j < n_b_sh; j += 1 ) 
          {
          for ( int ss = 0; ss < basisMatrix->Rows(); ss += 1 ) 
            {
            (*basisMatrix)(ss,i*n_b_sh+j) = (*this->m_BasisMatrix)(ss,(i+1)*n_b_sh+j) - (*m_Gn0)[i+1]/(*m_Gn0)[0] * (*this->m_BasisMatrix)(ss,j);
            }
          }
        }
      if (this->m_EstimationType==Self::L1_DL)
        {
        if (this->m_BasisCombinationMatrix->Size()==0)
          utl::ReadMatrix(utl::CreateExpandedPath(utl::LearnedSPFDictionary_SH8_RA4_K250), *this->m_BasisCombinationMatrix);
        // basisMatrix *= (*this->m_BasisCombinationMatrix);
        MatrixPointer tmpMat (new MatrixType());
        utl::ProductUtlMM(*basisMatrix, *this->m_BasisCombinationMatrix, *tmpMat);
        basisMatrix = tmpMat;
        }
      }
    if (this->GetDebug())
      utl::PrintUtlMatrix(*basisMatrix, "basisMatrix_InSolver");

    if (this->m_EstimationType==Self::LS)
      {
      this->m_L2Solver->SetA(basisMatrix);
      }
    else if (this->m_EstimationType==Self::L1_2 || this->m_EstimationType==Self::L1_DL)
      {
      if (this->m_L1SolverType==Self::FISTA_LS)
        {
        this->m_L1FISTASolver->SetA(basisMatrix);
        // this->m_L1FISTASolver->Print(std::cout<<"this->m_L1FISTASolver 0 = ");
        }
      else if (this->m_L1SolverType==Self::SPAMS)
        this->m_L1SpamsSolver->SetA(basisMatrix);
      else 
        utlException(true, "wrong m_L1SolverType");
      }
    }
  else
    {
    std::cout << "Use adaptive scale for each voxel!" << std::endl << std::flush;
    // calculate scale image
    typedef SPFScaleFromMeanDiffusivityImageFilter<ScalarImageType, ScalarImageType>  ScaleFromMDfilterType;
    typename ScaleFromMDfilterType::Pointer scaleFromMDfilter = ScaleFromMDfilterType::New();
    scaleFromMDfilter->SetMD0(this->m_MD0);
    scaleFromMDfilter->SetTau(this->m_SamplingSchemeQSpace->GetTau());
    scaleFromMDfilter->SetIsOriginalBasis(this->m_IsOriginalBasis);
    scaleFromMDfilter->SetInput(this->m_MDImage);
    scaleFromMDfilter->SetInPlace(false); // no inplace
    scaleFromMDfilter->Update();
    this->m_ScaleImage = scaleFromMDfilter->GetOutput();
    }

  this->ComputeRegularizationWeight();
  if (this->GetDebug())
    utl::PrintUtlVector(*this->m_RegularizationWeight, "m_RegularizationWeight");

  STDVectorPointer qVector = this->m_SamplingSchemeQSpace->GetRadiusVector();
  if (this->m_EstimationType==Self::LS)
    {
    if (this->m_LambdaSpherical>0 || this->m_LambdaRadial>0)
      {
      MatrixPointer lamMat (new MatrixType());
      *lamMat = this->m_RegularizationWeight->GetDiagonalMatrix();
      this->m_L2Solver->SetLambda(lamMat);
      }
    }
  else if (this->m_EstimationType==Self::L1_2)
    {
    if (this->m_L1SolverType==Self::FISTA_LS)
      {
      this->m_L1FISTASolver->SetwForInitialization(this->m_RegularizationWeight);
      this->m_L1FISTASolver->Setw(this->m_RegularizationWeight);
      // this->m_L1FISTASolver->Setw(*this->m_RegularizationWeight * qVector->size());
      }
    else if (this->m_L1SolverType==Self::SPAMS)
      this->m_L1SpamsSolver->Setw(this->m_RegularizationWeight);
    }
  else if (this->m_EstimationType==Self::L1_DL)
    {
    if (this->m_L1SolverType==Self::FISTA_LS)
      {
      // give a small regularization for the initialization in L2Solver
      // this->m_L1FISTASolver->SetwForInitialization(*this->m_RegularizationWeight/this->m_RegularizationWeight->two_norm()*1e-7);
      this->m_L1FISTASolver->SetwForInitialization(utl::ToVector<double> (*this->m_RegularizationWeight % qVector->size()) );
      // tune the weight based on the size of DWI samples
      this->m_L1FISTASolver->Setw(utl::ToVector<double>(*this->m_RegularizationWeight % qVector->size()) );
      }
    else if (this->m_L1SolverType==Self::SPAMS)
      this->m_L1SpamsSolver->Setw(utl::ToVector<double>(*this->m_RegularizationWeight % qVector->size()));
    }

  this->InitializeThreadedLibraries();

  // this->m_L2Solver->Initialize();
  // MatrixType ls = this->m_L2Solver->GetLS();
  // utl::PrintUtlMatrix(ls,"LS");
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId )
{
  itkShowPositionThreadedLogger(this->GetDebug());
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  // Pointers
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer outputPtr = this->GetOutput();

  // iterator for the output image      
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputRegionForThread );
  ImageRegionConstIteratorWithIndex<InputImageType> inputIt(inputPtr, outputRegionForThread );
  ImageRegionIteratorWithIndex<MaskImageType> maskIt;
  ImageRegionIteratorWithIndex<ScalarImageType> scaleIt;
  if (this->IsMaskUsed())
    maskIt = ImageRegionIteratorWithIndex<MaskImageType>(this->m_MaskImage, outputRegionForThread);
  if (!IsImageEmpty(this->m_ScaleImage))
    scaleIt = ImageRegionIteratorWithIndex<ScalarImageType>(this->m_ScaleImage, outputRegionForThread);
  
  InputImagePixelType inputPixel;
  // OutputImageIndexType outputIndex;
  OutputImagePixelType outputPixel, outputZero;
  
  unsigned int numberOfCoeffcients = outputPtr->GetNumberOfComponentsPerPixel();;
  outputPixel.SetSize(numberOfCoeffcients);
  outputZero.SetSize(numberOfCoeffcients), outputZero.Fill(0.0);
  unsigned int numberOfDWIs = inputPtr->GetNumberOfComponentsPerPixel();
  inputPixel.SetSize(numberOfDWIs);
  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n_b_ra = this->m_RadialRank + 1;
  int n_b =  n_b_ra * n_b_sh;
  
  VectorType dwiPixel(numberOfDWIs+this->m_BasisMatrixForB0->Rows()), dwiPixel_est(numberOfDWIs+this->m_BasisMatrixForB0->Rows()), dwiPixel_first(numberOfDWIs), coef(numberOfCoeffcients), coef_first;
  InputImageIndexType index;


  Pointer selfClone = this->Clone();
  // change selfClone->m_ThreadID and generate threadIDStr
  selfClone->m_ThreadID = threadId;
  std::string threadIDStr =  selfClone->ThreadIDToString();
  if (this->GetDebug())
    {
    std::ostringstream msg;
    selfClone->Print(msg << threadIDStr <<"selfClone = ");
    this->WriteLogger(msg.str());
    }
  

  STDVectorPointer qVector = selfClone->m_SamplingSchemeQSpace->GetRadiusVector();

  if (!this->IsAdaptiveScale())
    {
    selfClone->ComputeRadialVectorForE0InDWI();
    selfClone->ComputeRadialVectorForE0InBasis();
    }

  MatrixPointer basisMatrix(new MatrixType());
  if (!this->m_IsAnalyticalB0 && !this->IsAdaptiveScale())
    *basisMatrix = utl::ConnectUtlMatrix(*selfClone->m_BasisMatrix, *utl::ToMatrix<double>(*selfClone->m_BasisMatrixForB0 % selfClone->m_B0Weight), true);
  
  for (inputIt.GoToBegin(), outputIt.GoToBegin(), maskIt.GoToBegin(), scaleIt.GoToBegin();
    !inputIt.IsAtEnd(); 
    progress.CompletedPixel(), ++inputIt, ++outputIt, ++maskIt, ++scaleIt) 
    {

    if (this->IsMaskUsed() && maskIt.Get()<=1e-8)
      {
      outputIt.Set(outputZero);
      continue;
      }

    inputPixel=inputIt.Get();
    if (inputPixel.GetSquaredNorm()<=1e-8)
      {
      outputIt.Set(outputZero);
      continue;
      }

    index = inputIt.GetIndex();
    if (this->GetDebug())
      {
      std::ostringstream msg;
      msg << "\n" << threadIDStr << "index = " << index << std::endl << std::flush;
      this->WriteLogger(msg.str());
      }

    for ( int i = 0; i < numberOfDWIs; i += 1 ) 
      dwiPixel[i] = inputPixel[i];

    if (this->IsAdaptiveScale())
      {
      double scale = scaleIt.Get();
      if (scale<=0)
        {
        // use scaleImage as a mask
        outputIt.Set(outputZero);
        continue;
        }

      selfClone->SetBasisScale(scale);
      if (selfClone->m_BasisMatrix->Rows()==0)
        {
        selfClone->ComputeBasisMatrix();
        if (!selfClone->m_IsAnalyticalB0)
          selfClone->ComputeBasisMatrixForB0();
        else 
          {
          selfClone->ComputeRadialVectorForE0InDWI();
          selfClone->ComputeRadialVectorForE0InBasis();
          }
        }

      if (!this->m_IsAnalyticalB0)
        {
        basisMatrix=MatrixPointer( new MatrixType() );
        *basisMatrix = utl::ConnectUtlMatrix(*selfClone->m_BasisMatrix, *utl::ToMatrix<double>((*selfClone->m_BasisMatrixForB0)%this->m_B0Weight), true);
        }
      else
        {
        if (basisMatrix->Size()==0 || this->m_EstimationType==Self::L1_DL) // resize it for L1_DL
          basisMatrix=MatrixPointer( new MatrixType(selfClone->m_BasisMatrix->Rows(), n_b_sh*(n_b_ra-1)) );
        utlException((*selfClone->m_Gn0)[0]==0, "it should be not zero!, (*selfClone->m_Gn0)[0]="<< (*selfClone->m_Gn0)[0]);
        // utlPrintVar2 (selfClone->_is_b0_analytical, selfClone->m_RadialRank, R_N_0);


        // utl::Tic(std::cout<<"index start 1");
        double *selfBasisMatrix_data = selfClone->m_BasisMatrix->GetData();
        double *basisMatrix_data = basisMatrix->GetData();
        int index_B=0, index_selfB=0, index_selfB_0=0;
        for ( int ss = 0; ss < basisMatrix->Rows(); ss += 1 ) 
          {
          for ( int i = 0; i < this->m_RadialRank; i += 1 ) 
            {
            if (i==0)
              index_selfB += n_b_sh;
            else
              index_selfB_0 -= n_b_sh;
            for ( int j = 0; j < n_b_sh; j += 1 ) 
              {
              basisMatrix_data[index_B] = selfBasisMatrix_data[index_selfB] - (*selfClone->m_Gn0)[i+1]/(*selfClone->m_Gn0)[0] * selfBasisMatrix_data[index_selfB_0];
              // basisMatrix(ss,i*n_b_sh+j) = (*selfClone->m_BasisMatrix)(ss,(i+1)*n_b_sh+j) - (*selfClone->m_Gn0)[i+1]/(*selfClone->m_Gn0)[0] * (*selfClone->m_BasisMatrix)(ss,j);
              index_B++;
              index_selfB++;
              index_selfB_0++;
              }
            }
          index_selfB_0 += n_b_sh*this->m_RadialRank;
          }
        // utl::Toc();

        // utl::Tic(std::cout<<"() start 1");
        // for ( int i = 0; i < this->m_RadialRank; i += 1 ) 
        //   {
        //   double firstTerm = (*selfClone->m_Gn0)[i+1]/(*selfClone->m_Gn0)[0];
        //   for ( int j = 0; j < n_b_sh; j += 1 ) 
        //     {
        //     int index_B = i*n_b_sh+j;
        //     int index_selfB = (i+1)*n_b_sh+j;
        //     for ( int ss = 0; ss < basisMatrix->Rows(); ss += 1 ) 
        //       {
        //       (*basisMatrix)(ss,index_B) = (*selfClone->m_BasisMatrix)(ss,index_selfB) - firstTerm * (*selfClone->m_BasisMatrix)(ss,j);
        //       }
        //     }
        //   }
        // utl::Toc();

        }

      if (this->m_EstimationType==Self::L1_DL)
        {
        MatrixPointer tmpMat (new MatrixType());
        utl::ProductUtlMM(*basisMatrix, *this->m_BasisCombinationMatrix, *tmpMat);
        *basisMatrix = *tmpMat;
        // utl::MatrixCopy(*tmpMat, *basisMatrix, 1.0, 'N');
        // basisMatrix = tmpMat;
        // basisMatrix *= (*this->m_BasisCombinationMatrix);
        }

      if (this->GetDebug())
        {
        std::ostringstream msg;
        utl::PrintUtlMatrix(*basisMatrix, "basisMatrix_InSolver", " ", msg << threadIDStr);
        this->WriteLogger(msg.str());
        }

      if (this->m_EstimationType==Self::LS)
        {
        selfClone->m_L2Solver->SetA(basisMatrix);
        }
      else if (this->m_EstimationType==Self::L1_2 || this->m_EstimationType==Self::L1_DL)
        {
        if (this->m_L1SolverType==Self::FISTA_LS)
          selfClone->m_L1FISTASolver->SetA(basisMatrix);
        else if (this->m_L1SolverType==Self::SPAMS)
          selfClone->m_L1SpamsSolver->SetA(basisMatrix);
        }
      }

    if (this->m_IsAnalyticalB0)
      {
      for ( int i = 0; i < dwiPixel.Size(); i += 1 ) 
        dwiPixel_first[i] = dwiPixel[i] - (*selfClone->m_G0DWI)[i];

      if (this->m_EstimationType==Self::LS)
        {
        selfClone->m_L2Solver->Setb(VectorPointer(new VectorType(dwiPixel_first)));
        selfClone->m_L2Solver->Solve();
        coef_first = selfClone->m_L2Solver->Getx();
        // MatrixType ls = selfClone->m_L2Solver->GetLS();
        // utl::PrintUtlMatrix(ls, "LS");
        // VectorType coef_test = ls*dwiPixel_first;
        // utl::PrintUtlVector(coef_test, "coef_test");
        }
      else if (this->m_EstimationType==Self::L1_2 || this->m_EstimationType==Self::L1_DL)
        {
        if (this->m_L1SolverType==Self::FISTA_LS)
          {
          selfClone->m_L1FISTASolver->Setb(VectorPointer(new VectorType(dwiPixel_first)));
          // selfClone->m_L1FISTASolver->SetDebug(this->GetDebug());
          // selfClone->m_L1FISTASolver->Print(std::cout<<"selfClone->m_L1FISTASolver = ");
          selfClone->m_L1FISTASolver->Solve();
          coef_first = selfClone->m_L1FISTASolver->Getx();
          }
        else if (this->m_L1SolverType==Self::SPAMS)
          {
          selfClone->m_L1SpamsSolver->Setb(VectorPointer(new VectorType(dwiPixel_first)));
          selfClone->m_L1SpamsSolver->Solve();
          coef_first = selfClone->m_L1SpamsSolver->Getx();
          }

        if (this->GetDebug())
          {
          std::ostringstream msg;
          if (selfClone->m_L1SpamsSolver)
            {
            msg << threadIDStr << "use m_L1SpamsSolver" << std::endl << std::flush;
            selfClone->m_L1SpamsSolver->Print(msg << threadIDStr<<"this->m_L1QPSolver = ");
            double func = selfClone->m_L1SpamsSolver->EvaluateCostFunction();
            msg << threadIDStr << "func spams = " << func << std::endl << std::flush;
            }
          if (selfClone->m_L1FISTASolver)
            {
            msg << threadIDStr << "use m_L1FISTASolver" << std::endl << std::flush;
            selfClone->m_L1FISTASolver->Print(msg << threadIDStr<<"this->m_L1FISTASolver = ");
            std::vector<double> funcVec = selfClone->m_L1FISTASolver->GetCostFunction();
            utl::PrintVector(funcVec, "func FISTA", " ", msg << threadIDStr);
            }
          this->WriteLogger(msg.str());
          }
        }

      if (this->m_EstimationType==Self::L1_DL)
        {
        VectorType coef_first_tmp;
        utl::ProductUtlMv(*this->m_BasisCombinationMatrix, coef_first,coef_first_tmp);
        for ( int i = 0; i < coef_first_tmp.Size(); i += 1 ) 
          coef[i+n_b_sh] = coef_first_tmp[i];
        }
      else
        {
        for ( int i = 0; i < coef_first.Size(); i += 1 ) 
          coef[i+n_b_sh] = coef_first[i];
        }

      int jj=0;
      for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
        {
        for ( int m = -l; m <= l; m += 1 ) 
          {
          double sum_tmp = 0;
          for ( int nn = 1; nn <= this->m_RadialRank; nn += 1 ) 
            {
            int index_j = this->GetIndexJ(nn,l,m);
            sum_tmp += coef[index_j] * (*selfClone->m_Gn0)[nn];
            }
          if (l==0)
            coef[jj] = (std::sqrt(4*M_PI) - sum_tmp)/(*selfClone->m_Gn0)[0];
          else
            coef[jj] = -sum_tmp/(*selfClone->m_Gn0)[0];
          jj++;
          }
        }

      }
    else
      {
      for ( int i = 0; i < selfClone->m_BasisMatrixForB0->Rows(); i += 1 ) 
        dwiPixel[i+numberOfDWIs] = selfClone->m_B0Weight;
      // utl::PrintUtlVector(dwiPixel, "dwiPixel");
      // utl::PrintUtlMatrix(selfClone->m_BasisMatrixForB0, "selfClone->m_BasisMatrixForB0");

      if (this->m_EstimationType==Self::LS)
        {
        selfClone->m_L2Solver->Setb(VectorPointer(new VectorType(dwiPixel)));
        selfClone->m_L2Solver->Solve();
        coef = selfClone->m_L2Solver->Getx();
        }
      else if (this->m_EstimationType==Self::L1_2 || this->m_EstimationType==Self::L1_DL)
        {
        if (this->m_L1SolverType==Self::FISTA_LS)
          {
          selfClone->m_L1FISTASolver->Setb(VectorPointer(new VectorType(dwiPixel)));
          selfClone->m_L1FISTASolver->Solve();
          coef = selfClone->m_L1FISTASolver->Getx();
          }
        else if (this->m_L1SolverType==Self::SPAMS)
          {
          selfClone->m_L1SpamsSolver->Setb(VectorPointer(new VectorType(dwiPixel)));
          selfClone->m_L1SpamsSolver->Solve();
          coef = selfClone->m_L1SpamsSolver->Getx();
          }
        }

      }

    if (this->GetDebug())
      {
      std::ostringstream msg;
      utl::PrintUtlVector(dwiPixel, "dwiPixel", " ", msg << threadIDStr);
      utl::PrintUtlVector(coef, "coef", " ", msg << threadIDStr);
      if (this->m_IsAnalyticalB0)
        {
        utl::PrintUtlVector(dwiPixel_first, "dwiPixel_first", " ", msg << threadIDStr);
        utl::PrintUtlVector(coef_first, "coef_first", " ", msg << threadIDStr);
        utl::ProductUtlMv(*selfClone->m_BasisMatrix, coef, dwiPixel_est);
        }
      else
        {
        utl::ProductUtlMv(*basisMatrix, coef, dwiPixel_est);
        }
      utl::PrintUtlVector(dwiPixel_est, "dwiPixel_est", " ", msg << threadIDStr);
      VectorType diff = (dwiPixel - dwiPixel_est); diff.ElementAbsolute();
      utl::PrintUtlVector(diff, "dwiPixel_diff", " ", msg << threadIDStr);
      if (selfClone->m_BasisMatrixForB0->Size()==0)
        selfClone->ComputeBasisMatrixForB0();
      VectorType dwiPixel_b0_est;
      utl::ProductUtlMv(*selfClone->m_BasisMatrixForB0, coef, dwiPixel_b0_est);
      utl::PrintUtlVector(dwiPixel_b0_est, "dwiPixel_b0_est", " ", msg << threadIDStr);
      this->WriteLogger(msg.str());
      }

    for ( int i = 0; i < numberOfCoeffcients; i += 1 ) 
      outputPixel[i] = coef[i];

    outputIt.Set(outputPixel);
    }
}

template< class TInputImage, class TOutputImage >
void
SphericalPolarFourierImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

}


#endif 

