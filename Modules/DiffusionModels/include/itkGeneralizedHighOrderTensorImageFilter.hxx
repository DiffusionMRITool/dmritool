/**
 *       @file  itkGeneralizedHighOrderTensorImageFilter.hxx
 *      @brief  
 *
 *
 *     @author  Dr. Jian Cheng (JC), jian.cheng.1983@gmail.com
 *
 *   @internal
 *     Created  "04-22-2013
 *    Revision  1.0
 *    Compiler  gcc/g++
 *     Company  IDEA@UNC-CH
 *   Copyright  Copyright (c) 2013, Jian Cheng
 *
 * =====================================================================================
 */

#ifndef __itkGeneralizedHighOrderTensorImageFilter_hxx
#define __itkGeneralizedHighOrderTensorImageFilter_hxx

#include "itkGeneralizedHighOrderTensorImageFilter.h"
#include "itkProgressReporter.h"

namespace itk
{
template< class TInputImage, class TOutputImage >
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::GeneralizedHighOrderTensorImageFilter() : Superclass()
{
}

template< class TInputImage, class TOutputImage >
double
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::ComputeScale(const bool setScale)
{
  utlShowPosition(this->GetDebug());
  double scale = -1;
  STDVectorPointer bVector = this->m_SamplingSchemeQSpace->GetBVector();
  MatrixPointer qOrientations = this->m_SamplingSchemeQSpace->GetOrientationsSpherical();
  // utl::PrintVector(*this->m_SamplingSchemeQSpace->GetBVector(), "this->m_SamplingSchemeQSpace->GetBVector()");
  // utl::PrintUtlMatrix(*qOrientations, "qOrientations");
  utlGlobalException(bVector->size()==0, "no b values");
  utlGlobalException(qOrientations->Rows()==0, "no gradients");
  double bMax = *std::max_element(this->m_SamplingSchemeQSpace->GetBVector()->begin(), this->m_SamplingSchemeQSpace->GetBVector()->end());
  scale = 0.5*bMax/(4*M_PI*M_PI*this->m_SamplingSchemeQSpace->GetTau());
  if (setScale)
    {
    this->SetBasisScale(scale);
    }
  if (this->GetDebug())
    std::cout << "m_BasisScale = " << this->m_BasisScale << std::endl;
  return scale;

}

template< class TInputImage, class TOutputImage >
std::vector<int>
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::DimToRank ( const int dimm ) const
{
  std::vector<int> result;
  int radialRank=-1, shRank=-1;
  for ( int radialRank = 1; radialRank <= 10; radialRank += 1 ) 
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
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::RankToDim (const bool is_radial, const int radialRank, const int shRank) const
{
  int radialRank_real = radialRank>=0?radialRank:this->m_RadialRank;
  int shRank_real = shRank>=0?shRank:this->m_SHRank;
  if (is_radial)
    return radialRank_real;
  else
    return (shRank_real + 1)*(shRank_real + 2)/2*(radialRank_real);
}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::VerifyInputParameters() const
{
  utlShowPosition(this->GetDebug());
  Superclass::VerifyInputParameters();
  utlGlobalException(this->m_RadialRank<=0, "m_RadialRank should be no less than 1");
  utlGlobalException(this->m_EstimationType==Superclass::L1_2, "TODO");
  utlGlobalException(this->m_EstimationType==Superclass::L1_DL, "TODO");
}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::ComputeRadialMatrix ()
{
  utlShowPosition(this->GetDebug());

  int n_s, n_b;
  n_s = this->m_SamplingSchemeQSpace->GetBVector()->size();
  utlGlobalException( n_s==0, "no b vector");

  if (this->GetDebug())
    std::cout << "m_BasisScale = " << this->m_BasisScale << std::endl;

  const STDVectorPointer bVector = this->m_SamplingSchemeQSpace->GetBVector();

  // NOTE: b = 4\pi^2 * _tau * q^2
  STDVectorType qVector(*bVector);
  for ( int i = 0; i < bVector->size(); i += 1 ) 
    {
    qVector[i] = std::sqrt((*bVector)[i]/(4*M_PI*M_PI*this->m_SamplingSchemeQSpace->GetTau())); 
    }


  // the basis of order N has N+1 terms (1,...,N), NOTE: does not start from 0
  n_b = this->m_RadialRank;


  if(this->GetDebug())
    std::cout << "Generating the "<< n_s << "x" << n_b << " RadialMatrix...\n";

  MatrixPointer B(new MatrixType(n_s,n_b));

  for ( int js = 0; js < n_s; js += 1 ) 
    {
    double x_temp = qVector[js] / std::sqrt(this->m_BasisScale);
    for ( int ib = 0; ib < n_b; ib += 1 ) 
      {
      (*B)(js,ib) = std::pow(x_temp, 2.0*(ib+1));
      }
    }


  if(this->GetDebug()) 
    {
    std::cout << "Generated the "<< n_s << "x" << n_b << " RadialMatrix...\n";
    utl::PrintUtlMatrix(*B,"RadialMatrix");
    }

  this->m_BasisRadialMatrix = B;

}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::ComputeBasisMatrix() 
{
  utlShowPosition(this->GetDebug());

  if (this->m_BasisSHMatrix->Rows()==0)
    this->ComputeSHMatrix();
  MatrixPointer B_sh = this->m_BasisSHMatrix;

  if (this->m_BasisRadialMatrix->Rows()==0)
    this->ComputeRadialMatrix();
  MatrixPointer B_ra = this->m_BasisRadialMatrix;
  
  MatrixPointer qOrientations = this->m_SamplingSchemeQSpace->GetOrientationsSpherical();

  int bVector_size = this->m_SamplingSchemeQSpace->GetBVector()->size();
  int grad_size = qOrientations->Rows();

  if (this->GetDebug())
    {
    std::cout << "this->m_SamplingSchemeQSpace->GetBVector()->size() = " << this->m_SamplingSchemeQSpace->GetBVector()->size() << std::endl;
    std::cout << "qOrientations->Rows() = " << qOrientations->Rows() << std::endl;
    std::cout << "bVector_size = " << bVector_size << std::endl;
    std::cout << "grad_size = " << grad_size << std::endl;
    }
  utlException(bVector_size!=grad_size, "bVector_size and grad_size should keep the same size");
  int n_s = bVector_size;

  // the basis has two parts 
  int n_b_ra = this->m_RadialRank;
  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n_b =  n_b_ra * n_b_sh;

  if(this->GetDebug())
    {
    std::cout << "n_b_ra = " << n_b_ra  << std::endl;
    std::cout << "n_b_sh = " << n_b_sh  << std::endl;
    std::cout << "n_b = " << n_b  << std::endl;
    std::cout << "Generating the "<< n_s << "x" << n_b << " basis matrix...\n";
    }

  MatrixPointer B(new MatrixType(n_s, n_b));
  utlException(B_sh->Columns()!=n_b_sh, "the SHMatrix does not have the right width. B_sh=" << *B_sh << ", n_b_sh="<< n_b_sh);
  utlException(B_sh->Rows()!=B_ra->Rows(), "the SHMatrix and the RadialMatrix do not have the same samples. B_sh="<< *B_sh << ", B_ra=" << *B_ra);
  utlException(B_sh->Rows()!=B->Rows(), "the SHMatrix and the basisMatrix do not have the same samples. B_sh="<< *B_sh << ", B="<< *B);

  utlException(B_ra->Columns()!=n_b_ra, "the RadialMatrix does not have the right width");
  for ( int i = 0; i < n_b_ra; i += 1 ) 
    for ( int j = 0; j < n_b_sh; j += 1 ) 
      for ( int k = 0; k < n_s; k += 1 ) 
        {
        (*B)(k,i*n_b_sh+j) = (*B_ra)(k,i) * (*B_sh)(k,j); 
        }

  if(this->GetDebug()) 
    {
    utl::PrintUtlMatrix(*B,"BasisMatrix");
    }

  this->m_BasisMatrix = B;
}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::ComputeRegularizationWeight( )
{
  // utlShowPosition(this->GetDebug());

  int n_b_sh = (this->m_SHRank+1)*(this->m_SHRank+2)/2;
  int n_b_ra = this->m_RadialRank;

  this->m_RegularizationWeight = VectorPointer(new VectorType(this->RankToDim()));
  for ( int i = 0; i <= this->m_RadialRank-1; i += 1 ) 
    {
    int j = 0;
    for ( int l = 0; l <= this->m_SHRank; l += 2 ) 
      {
      for ( int m = -l; m <= l; m += 1 ) 
        {
        (*this->m_RegularizationWeight)(i*n_b_sh+j) = this->m_LambdaSpherical*l*l*(l+1)*(l+1) + this->m_LambdaRadial*i*i*(i+1)*(i+1);
        j++;
        }
      }
    }
}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  utlShowPosition(this->GetDebug());
  ComputeScale();
  ComputeBasisMatrix();
  this->VerifyInputParameters();
  this->m_L2Solver = L2SolverType::New();
  this->m_L2Solver->SetA(this->m_BasisMatrix);
  if (this->m_LambdaSpherical>0 || this->m_LambdaRadial>0)
    {
    this->ComputeRegularizationWeight();
    MatrixPointer mat(new MatrixType(this->m_RegularizationWeight->Size(), this->m_RegularizationWeight->Size()));
    mat->SetDiagonal(*this->m_RegularizationWeight);
    this->m_L2Solver->SetLambda(mat);
    }
  // if (this->GetDebug())
  //   {
  //   this->m_L2Solver->Initialize();
  //   MatrixType ls = this->m_L2Solver->GetLS();
  //   utl::PrintUtlMatrix(ls,"LS");
  //   }
}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,ThreadIdType threadId )
{
  utlShowPosition(this->GetDebug());
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  // Pointers
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer outputPtr = this->GetOutput();
  
  Pointer selfClone = this->Clone();

  // iterator for the output image      
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(outputPtr, outputRegionForThread );
  ImageRegionConstIteratorWithIndex<InputImageType> inputIt(inputPtr, outputRegionForThread );
  ImageRegionIteratorWithIndex<MaskImageType> maskIt;
  if (this->IsMaskUsed())
    maskIt = ImageRegionIteratorWithIndex<MaskImageType>(this->m_MaskImage, outputRegionForThread);
  
  InputImagePixelType inputPixel;
  // OutputImageIndexType outputIndex;
  OutputImagePixelType outputPixel;
  
  unsigned int numberOfCoeffcients = outputPtr->GetNumberOfComponentsPerPixel();;
  outputPixel.SetSize(numberOfCoeffcients);
  unsigned int numberofDWIs = inputPtr->GetNumberOfComponentsPerPixel();
  inputPixel.SetSize(numberofDWIs);
  
  inputIt.GoToBegin();
  outputIt.GoToBegin();
  VectorType dwiPixel(numberofDWIs), coef(numberOfCoeffcients);
  while( !inputIt.IsAtEnd() )
    {
    if (!this->IsMaskUsed() || (this->IsMaskUsed() && maskIt.Get()>0))
      {
      inputPixel=inputIt.Get();
      for ( int i = 0; i < numberofDWIs; i += 1 ) 
        dwiPixel[i] = -std::log(inputPixel[i]);
      selfClone->m_L2Solver->Setb(VectorPointer(new VectorType(dwiPixel)));
      // outputIndex=outputIt.GetIndex();
      // std::cout << "index="<<outputIndex << std::endl << std::flush;
      selfClone->m_L2Solver->Solve();
      coef = selfClone->m_L2Solver->Getx();
      for ( int i = 0; i < numberOfCoeffcients; i += 1 ) 
        outputPixel[i] = coef[i];
      // utl::PrintContainer(inputPixel.GetDataPointer(), inputPixel.GetDataPointer()+inputPixel.GetSize(), "dwi");
      // utl::PrintContainer(outputPixel.GetDataPointer(), inputPixel.GetDataPointer()+outputPixel.GetSize(), "coef");
      }
    else
      outputPixel.Fill(0.0);

    outputIt.Set(outputPixel);
    progress.CompletedPixel();    

    if (this->IsMaskUsed())
      ++maskIt;
    ++outputIt;
    ++inputIt;
    }
}

template< class TInputImage, class TOutputImage >
void
GeneralizedHighOrderTensorImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

}

#endif 
