#ifndef __itkSlowPolyLineParametricPath_hxx
#define __itkSlowPolyLineParametricPath_hxx

#include "itkSlowPolyLineParametricPath.h"
#include <math.h>



namespace itk
{

//template<unsigned int VDimension>
//typename SlowPolyLineParametricPath<VDimension>::VectorType
//SlowPolyLineParametricPath<VDimension>
//::EvaluateDerivative(const InputType & input) const
//{
//}



/**
 * Constructor
 */
template <unsigned int VDimension>
SlowPolyLineParametricPath<VDimension>
::SlowPolyLineParametricPath()
{
  this->SetDefaultInputStepSize( 0.3 );
  m_UseCentralDifference=true;
}

template< unsigned int VDimension >
typename PolyLineParametricPath< VDimension >::OutputType
SlowPolyLineParametricPath< VDimension >
::Evaluate(const InputType & input) const
{
  // Handle the endpoint carefully, since there is no following vertex
  const VertexListType* vertexList = this->GetVertexList();
  const InputType endPoint = static_cast< InputType >( vertexList->Size() - 1 );
  if ( input > endPoint || itk::Math::FloatAlmostEqual( input, endPoint ) )
    {
    return vertexList->ElementAt(vertexList->Size() - 1); // the last vertex
    }
  if (input<0 || itk::Math::FloatAlmostEqual( input, 0.0 ))
    {
    return vertexList->ElementAt(0); // the first vertex
    }

  const VertexType vertex0 = vertexList->ElementAt( int(input)<0.0 ? 0 : int(input)  );
  const VertexType vertex1 = vertexList->ElementAt( int(input)+1>=endPoint ? endPoint-1 : int(input)+1);

  const double fractionOfLineSegment = input - int(input);

  const PointType outputPoint = vertex0 + ( vertex1 - vertex0 ) * fractionOfLineSegment;

  // For some stupid reason, there is no easy way to cast from a point to a
  // continuous index.
  OutputType output;
  for ( unsigned int i = 0; i < VDimension; i++ )
    {
    output[i] = outputPoint[i];
    }

  return output;
}

template<unsigned int VDimension>
typename PolyLineParametricPath<VDimension>::VectorType
SlowPolyLineParametricPath<VDimension>
::EvaluateDerivative(const InputType & input, bool isDerivativeNormalizedByDistance) const
{
  if (m_UseCentralDifference)
    {
    //Get next integral time-point
    const InputType nextTimepoint = std::floor(input + 1.0);

    //Get previous integral time-point
    const InputType previousTimepoint = std::floor(input - 1.0);

    //Calculate the continuous index for both points
    const ContinuousIndexType nextIndex = this->Evaluate(nextTimepoint);
    const ContinuousIndexType prevIndex = this->Evaluate(previousTimepoint);

    //For some reason, there's no way to convert ContinuousIndexType to VectorType
    VectorType partialDerivatives;
    for (unsigned int i = 0; i < VDimension; ++i)
      {
      partialDerivatives[i] = (nextIndex[i] - prevIndex[i]);
      }

    if (isDerivativeNormalizedByDistance)
      {
      double dist = nextIndex.EuclideanDistanceTo(prevIndex);
      if (dist>1e-10)
        partialDerivatives /= dist;
      }

    return partialDerivatives;
    }
  else
    return Superclass::EvaluateDerivative(input);
}

template<unsigned int VDimension>
typename SlowPolyLineParametricPath<VDimension>::OffsetType
SlowPolyLineParametricPath<VDimension>
::IncrementInput(InputType & input) const
{
  int         iterationCount;
  bool        tooSmall;
  bool        tooBig;
  InputType   inputStepSize;
  InputType   finalInputValue;
  OffsetType  offset;
  IndexType   currentImageIndex;
  IndexType   nextImageIndex;
  IndexType   finalImageIndex;
  
  iterationCount    = 0;
  inputStepSize     = this->GetDefaultInputStepSize();

  // Are we already at (or past) the end of the input?
  finalInputValue   = this->EndOfInput();
  currentImageIndex = this->EvaluateToIndex( input );
  finalImageIndex   = this->EvaluateToIndex( finalInputValue );
  offset            = finalImageIndex - currentImageIndex;
  if(  ( offset == this->GetZeroOffset() && input != this->StartOfInput() )  ||
       ( input >=finalInputValue )  )
    {
    return this->GetZeroOffset();
    }
  
  do
    {
    if( iterationCount++ > 10000 ) {itkExceptionMacro(<<"Too many iterations");}
    
    nextImageIndex    = this->EvaluateToIndex( input + inputStepSize );
    offset            = nextImageIndex - currentImageIndex;
    
    tooBig = false;
    tooSmall = ( offset == this->GetZeroOffset() );
    if( tooSmall )
      {
      // increase the input step size, but don't go past the end of the input
      inputStepSize *= 2;
      if(  (input + inputStepSize) >= finalInputValue  ){
        //inputStepSize = finalInputValue - input;
        inputStepSize += this->GetDefaultInputStepSize();
      }
    }
    else
      {
      // Search for an offset dimension that is too big
      for( unsigned int i=0; i<VDimension && !tooBig; i++ )
        {
        tooBig = ( offset[i] >= 2 || offset[i] <= -2 );
        }
      
      if( tooBig ){
        //inputStepSize /= 1.5;
        inputStepSize -= (this->GetDefaultInputStepSize()/0.5);
      }
    }
  }
  while( tooSmall || tooBig );
  
  input += inputStepSize;
  return offset;
}


/**
 * Standard "PrintSelf" method
 */
template <unsigned int VDimension>
void
SlowPolyLineParametricPath<VDimension>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}



} // end namespaceitk

#endif
