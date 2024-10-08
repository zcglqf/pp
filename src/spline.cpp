// references:
//  https://github.com/InsightSoftwareConsortium/ITK/blob/master/Examples/RegistrationITKv4/DeformableRegistration12.cxx


#include "itkImageRegistrationMethodv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSBOptimizerv4.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageFilter.h"


#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#define DEBUG true

// Test data: 
//     fixed image: RARE64.nii
//     moving image: AtlasHead64.nii
//     output name : whatevername.nii
//     command: argv[0]  fullPathOfRARE64.nii fullPathOfAtlasHead64.nii fullPathOfwhatevername.nii
//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
    using Self = CommandIterationUpdate;
    using Superclass = itk::Command;
    using Pointer = itk::SmartPointer<Self>;
    itkNewMacro(Self);

protected:
    CommandIterationUpdate() = default;

public:
    using OptimizerType = itk::LBFGSBOptimizerv4;
    using OptimizerPointer = const OptimizerType*;

    void
        Execute(itk::Object* caller, const itk::EventObject& event) override
    {
        Execute((const itk::Object*)caller, event);
    }

    void
        Execute(const itk::Object* object, const itk::EventObject& event) override
    {
        auto optimizer = static_cast<OptimizerPointer>(object);
        if (!(itk::IterationEvent().CheckEvent(&event)))
        {
            return;
        }
        std::cout << optimizer->GetCurrentIteration() << "   ";
        std::cout << optimizer->GetCurrentMetricValue() << "   ";
        std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
    }
};



int main(int argc, char* argv[])
{
    std::cout << "Hellow post processing!";
    if (argc < 4)
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " fixedImageFile  movingImageFile outputImagefile  ";
        std::cerr << " [differenceOutputfile] [differenceBeforeRegistration] ";
        std::cerr << " [deformationField] ";
        std::cerr << " [filenameForFinalTransformParameters] ";
        std::cerr << " [numberOfGridNodesInOneDimension] ";
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Press any key to continue:";
    std::cin.get();

    constexpr unsigned int ImageDimension = 2;
    using PixelType = float;

    using ImageType = itk::Image<PixelType, ImageDimension>;

    // read fixed [or target] image (I suppose mri image is the fixed )
    using FixedImageReaderType = itk::ImageFileReader<ImageType>;
    auto fixedImageReader = FixedImageReaderType::New();
    fixedImageReader->SetFileName(argv[1]);
    fixedImageReader->Update();
    ImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
    ImageType::RegionType   fixedRegion = fixedImage->GetBufferedRegion();
   
    // read moving [or reference] image
    using MovingImageReaderType = itk::ImageFileReader<ImageType>;
    auto movingImageReader = MovingImageReaderType::New();
    movingImageReader->SetFileName(argv[2]);
    movingImageReader->Update();
    ImageType::ConstPointer movingImage = movingImageReader->GetOutput();

    ///////////// Debug only
    if (DEBUG)
	{		
		/////////// Debug only
		using MinMaxFilterType = itk::MinimumMaximumImageFilter<ImageType>;
		MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
		minMaxFilter->SetInput(fixedImageReader->GetOutput());
		minMaxFilter->Update();

		PixelType minimumValue = minMaxFilter->GetMinimum();
		PixelType maximumValue = minMaxFilter->GetMaximum();

		// Step 2: Extract image properties
		ImageType::PointType originFixed = fixedImage->GetOrigin();
		ImageType::SpacingType spacingFixed = fixedImage->GetSpacing();
		ImageType::SizeType sizeFixed = fixedImage->GetLargestPossibleRegion().GetSize();

		// Step 3: Compute the coverage of the image plane in real-world coordinates
		ImageType::PointType coverageStartFixed = originFixed;
		ImageType::PointType coverageEndFixed;

		for (unsigned int i = 0; i < ImageDimension; ++i)
		{
			coverageEndFixed[i] = originFixed[i] + spacingFixed[i] * (sizeFixed[i] - 1);
		}

		// Output the results
		std::cout << "Image Origin: " << originFixed << std::endl;
		std::cout << "Image Spacing: " << spacingFixed << std::endl;
		std::cout << "Image Size: " << sizeFixed << std::endl;
		std::cout << "Image Coverage in Real-World Coordinates:" << std::endl;
		std::cout << "  Start: " << coverageStartFixed << std::endl;
		std::cout << "  End:   " << coverageEndFixed << std::endl;
		std::cout << "  Image size: " << coverageEndFixed - coverageStartFixed << std::endl;
		std::cout << "  Minimum Value: " << minimumValue << std::endl;
		std::cout << "  Maximum Value: " << maximumValue << std::endl;

		//using MinMaxFilterType = itk::MinimumMaximumImageFilter<ImageType>;
		//MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
		minMaxFilter->SetInput(movingImage);
		minMaxFilter->Update();

		minimumValue = minMaxFilter->GetMinimum();
		maximumValue = minMaxFilter->GetMaximum();

		// Step 2: Extract image properties
		ImageType::PointType originMoving = movingImage->GetOrigin();
		ImageType::SpacingType spacingMoving = movingImage->GetSpacing();
		ImageType::SizeType sizeMoving = movingImage->GetLargestPossibleRegion().GetSize();

		// Step 3: Compute the coverage of the image plane in real-world coordinates
		ImageType::PointType coverageStartMoving = originMoving;
		ImageType::PointType coverageEndMoving;

		for (unsigned int i = 0; i < ImageDimension; ++i)
		{
			coverageEndMoving[i] = originMoving[i] + spacingMoving[i] * (sizeMoving[i] - 1);
		}

		// Output the results
		std::cout << "Image Origin: " << originMoving << std::endl;
		std::cout << "Image Spacing: " << spacingMoving << std::endl;
		std::cout << "Image Size: " << sizeMoving << std::endl;
		std::cout << "Image Coverage in Real-World Coordinates:" << std::endl;
		std::cout << "  Start: " << coverageStartMoving << std::endl;
		std::cout << "  End:   " << coverageEndMoving << std::endl;
		std::cout << "  Image size: " << coverageEndMoving - coverageStartMoving << std::endl;
		std::cout << "  Minimum Value: " << minimumValue << std::endl;
		std::cout << "  Maximum Value: " << maximumValue << std::endl;

		/////////////// Debug only
	}

    ///////////////////////////////////////////////////////////////////////////////////
    // Gauss smoothing
    ///////////////////////////////////////////////////////////////////////////////////
    using SmoothingFilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
    SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();
    smoothingFilter->SetInput(fixedImage);

    // Set the standard deviation (sigma) for the Gaussian kernel
    smoothingFilter->SetSigma(0.5);  // Adjust the sigma value as needed
    smoothingFilter->Update();
    // using WriterType = itk::ImageFileWriter<OutputImageType>;
    /// Gasss smoothing

    const unsigned int     SpaceDimension = ImageDimension;
    constexpr unsigned int SplineOrder = 1;
    using CoordinateRepType = double;

    using TransformType =
        itk::BSplineTransform<CoordinateRepType, SpaceDimension, SplineOrder>;
    using RegistrationType =
        itk::ImageRegistrationMethodv4<ImageType, ImageType>;
    auto registration = RegistrationType::New();
    auto transform = TransformType::New();
    unsigned int numberOfGridNodesInOneDimension = 7;


    if (argc > 8)
    {
        numberOfGridNodesInOneDimension = std::stoi(argv[8]);
    }

    // Software Guide : BeginCodeSnippet
    TransformType::PhysicalDimensionsType fixedPhysicalDimensions;
    TransformType::MeshSizeType           meshSize;
    TransformType::OriginType             fixedOrigin;

    for (unsigned int i = 0; i < SpaceDimension; ++i)
    {
        fixedOrigin[i] = fixedImage->GetOrigin()[i];
        fixedPhysicalDimensions[i] =
            fixedImage->GetSpacing()[i] *
            static_cast<double>(
                fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1);
    }
    std::cout << "Fixed image spatial diensions: " << fixedPhysicalDimensions;
    std::cout << std::endl;

    meshSize.Fill(numberOfGridNodesInOneDimension - SplineOrder);

    transform->SetTransformDomainOrigin(fixedOrigin);
    transform->SetTransformDomainPhysicalDimensions(fixedPhysicalDimensions);
    transform->SetTransformDomainMeshSize(meshSize);
    transform->SetTransformDomainDirection(fixedImage->GetDirection());

    registration->SetInitialTransform(transform);
    registration->InPlaceOn();

    using ParametersType = TransformType::ParametersType;

    const unsigned int numberOfParameters = transform->GetNumberOfParameters();

    ParametersType parameters(numberOfParameters);

    parameters.Fill(0.0);

    transform->SetParameters(parameters);
    //  Software Guide : EndCodeSnippet

    using MetricType =
        itk::MattesMutualInformationImageToImageMetricv4<ImageType,
        ImageType>;
    auto metric = MetricType::New();
    metric->SetNumberOfHistogramBins(32);
    metric->SetUseMovingImageGradientFilter(false);
    metric->SetUseFixedImageGradientFilter(false);
    metric->SetUseSampledPointSet(false);

    using OptimizerType = itk::LBFGSBOptimizerv4;
    auto optimizer = OptimizerType::New();

    // Software Guide : BeginCodeSnippet
    const unsigned int numParameters = transform->GetNumberOfParameters();
    OptimizerType::BoundSelectionType boundSelect(numParameters);
    OptimizerType::BoundValueType     upperBound(numParameters);
    OptimizerType::BoundValueType     lowerBound(numParameters);

    boundSelect.Fill(0);
    upperBound.Fill(0.0);
    lowerBound.Fill(0.0);

    optimizer->SetBoundSelection(boundSelect);
    optimizer->SetUpperBound(upperBound);
    optimizer->SetLowerBound(lowerBound);

    optimizer->SetCostFunctionConvergenceFactor(1.e7);
    optimizer->SetGradientConvergenceTolerance(1e-35);
    optimizer->SetNumberOfIterations(20);
    optimizer->SetMaximumNumberOfFunctionEvaluations(20);
    optimizer->SetMaximumNumberOfCorrections(3);
    // Software Guide : EndCodeSnippet

    // Create the Command observer and register it with the optimizer.
    //
    auto observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);

    // One level registration is performed using the shrink factor 1 and
    // smoothing sigma 1
    //
    constexpr unsigned int numberOfLevels = 1;

    RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
    shrinkFactorsPerLevel.SetSize(1);
    shrinkFactorsPerLevel[0] = 1;

    RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
    smoothingSigmasPerLevel.SetSize(1);
    smoothingSigmasPerLevel[0] = 0;

    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetNumberOfLevels(numberOfLevels);
    registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
    registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);


    // Add time and memory probes
    itk::TimeProbesCollectorBase   chronometer;
    itk::MemoryProbesCollectorBase memorymeter;

    std::cout << std::endl << "Starting Registration" << std::endl;

    try
    {
        memorymeter.Start("Registration");
        chronometer.Start("Registration");

        registration->Update();

        chronometer.Stop("Registration");
        memorymeter.Stop("Registration");

        const OptimizerType::ConstPointer outputOptimizer =
            dynamic_cast<const OptimizerType*>(registration->GetOptimizer());
        if (outputOptimizer.IsNotNull())
        {
            std::cout << "Optimizer stop condition = "
                << outputOptimizer->GetStopConditionDescription()
                << std::endl;
        }
        else
        {
            std::cerr << "Output optimizer is null." << std::endl;
            return EXIT_FAILURE;
        }
    }
    catch (const itk::ExceptionObject& err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    // While the registration filter is run, it updates the output transform
    // parameters with the final registration parameters
    OptimizerType::ParametersType finalParameters = transform->GetParameters();

    // Report the time and memory taken by the registration
    chronometer.Report(std::cout);
    memorymeter.Report(std::cout);

    using ResampleFilterType =
        itk::ResampleImageFilter<ImageType, ImageType>;

    auto resample = ResampleFilterType::New();

    resample->SetTransform(transform);
    resample->SetInput(movingImageReader->GetOutput());

    resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
    resample->SetOutputOrigin(fixedImage->GetOrigin());
    resample->SetOutputSpacing(fixedImage->GetSpacing());
    resample->SetOutputDirection(fixedImage->GetDirection());

    // This value is set to zero in order to make easier to perform
    // regression testing in this example. However, for didactic
    // exercise it will be better to set it to a medium gray value
    // such as 100 or 128.
    resample->SetDefaultPixelValue(0);

    using OutputPixelType = float;

    using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;

    using CastFilterType =
        itk::CastImageFilter<ImageType, OutputImageType>;

    using WriterType = itk::ImageFileWriter<OutputImageType>;


    auto writer = WriterType::New();
    auto caster = CastFilterType::New();


    writer->SetFileName(argv[3]);


    caster->SetInput(resample->GetOutput());
    writer->SetInput(caster->GetOutput());


    try
    {
        writer->Update();
    }
    catch (const itk::ExceptionObject& err)
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    using DifferenceFilterType =
        itk::SquaredDifferenceImageFilter<ImageType,
        ImageType,
        OutputImageType>;

    auto difference = DifferenceFilterType::New();

    auto writer2 = WriterType::New();
    writer2->SetInput(difference->GetOutput());


    // Compute the difference image between the
    // fixed and resampled moving image.
    if (argc > 4)
    {
        difference->SetInput1(fixedImageReader->GetOutput());
        difference->SetInput2(resample->GetOutput());
        writer2->SetFileName(argv[4]);
        try
        {
            writer2->Update();
        }
        catch (const itk::ExceptionObject& err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
    }


    // Compute the difference image between the
    // fixed and moving image before registration.
    if (argc > 5)
    {
        writer2->SetFileName(argv[5]);
        difference->SetInput1(fixedImageReader->GetOutput());
        difference->SetInput2(movingImageReader->GetOutput());
        try
        {
            writer2->Update();
        }
        catch (const itk::ExceptionObject& err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Generate the explicit deformation field resulting from
    // the registration.
    if (argc > 6)
    {

        using VectorType = itk::Vector<float, ImageDimension>;
        using DisplacementFieldType = itk::Image<VectorType, ImageDimension>;

        auto field = DisplacementFieldType::New();
        field->SetRegions(fixedRegion);
        field->SetOrigin(fixedImage->GetOrigin());
        field->SetSpacing(fixedImage->GetSpacing());
        field->SetDirection(fixedImage->GetDirection());
        field->Allocate();

        using FieldIterator = itk::ImageRegionIterator<DisplacementFieldType>;
        FieldIterator fi(field, fixedRegion);

        fi.GoToBegin();

        TransformType::InputPointType    fixedPoint;
        TransformType::OutputPointType   movingPoint;
        DisplacementFieldType::IndexType index;

        VectorType displacement;

        while (!fi.IsAtEnd())
        {
            index = fi.GetIndex();
            field->TransformIndexToPhysicalPoint(index, fixedPoint);
            movingPoint = transform->TransformPoint(fixedPoint);
            displacement = movingPoint - fixedPoint;
            fi.Set(displacement);
            ++fi;
        }

        using FieldWriterType = itk::ImageFileWriter<DisplacementFieldType>;
        auto fieldWriter = FieldWriterType::New();

        fieldWriter->SetInput(field);

        fieldWriter->SetFileName(argv[6]);
        try
        {
            fieldWriter->Update();
        }
        catch (const itk::ExceptionObject& excp)
        {
            std::cerr << "Exception thrown " << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Optionally, save the transform parameters in a file
    if (argc > 7)
    {
        std::ofstream parametersFile;
        parametersFile.open(argv[7]);
        parametersFile << finalParameters << std::endl;
        parametersFile.close();
    }

    return EXIT_SUCCESS;
}
