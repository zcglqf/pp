#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAffineTransform.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

constexpr unsigned int Dimension = 3;
using PixelType = double;
using ImageType = itk::Image<PixelType, Dimension>;

using TransformType = itk::AffineTransform<double, Dimension>;
using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
using MetricType = itk::MeanSquaresImageToImageMetricv4<ImageType, ImageType>;
using RegistrationType = itk::ImageRegistrationMethodv4<ImageType, ImageType, TransformType>;

using ReaderType = itk::ImageFileReader<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;

int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cerr << "Usage: " << argv[0]
            << " <FixedImage> <MovingImage> <MovingLabel> <OutputTransformedLabel>"
            << std::endl;
        return EXIT_FAILURE;
    }

    const char* fixedImageFileName = argv[1];
    const char* movingImageFileName = argv[2];
    const char* movingLabelFileName = argv[3];
    const char* outputTransformedLabelFileName = argv[4];

    // Read the fixed and moving images
    ReaderType::Pointer fixedImageReader = ReaderType::New();
    fixedImageReader->SetFileName(fixedImageFileName);
    fixedImageReader->Update();

    ReaderType::Pointer movingImageReader = ReaderType::New();
    movingImageReader->SetFileName(movingImageFileName);
    movingImageReader->Update();

    ReaderType::Pointer movingLabelReader = ReaderType::New();
    movingLabelReader->SetFileName(movingLabelFileName);
    movingLabelReader->Update();

    // Set up the registration
    RegistrationType::Pointer registration = RegistrationType::New();
    MetricType::Pointer metric = MetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    TransformType::Pointer transform = TransformType::New();

    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetFixedImage(fixedImageReader->GetOutput());
    registration->SetMovingImage(movingImageReader->GetOutput());
    registration->SetInitialTransform(transform);
    registration->InPlaceOn();

    // Configure the optimizer
    optimizer->SetLearningRate(4.0);
    optimizer->SetMinimumStepLength(0.001);
    optimizer->SetNumberOfIterations(200);

    // Perform the registration
    try
    {
        registration->Update();
    }
    catch (itk::ExceptionObject& err)
    {
        std::cerr << "ExceptionObject caught: " << err << std::endl;
        return EXIT_FAILURE;
    }

    // Get the final transform
    TransformType::ConstPointer finalTransform = registration->GetTransform();

    // Resample the moving label image
    using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(movingLabelReader->GetOutput());
    resampler->SetTransform(finalTransform);
    resampler->SetUseReferenceImage(true);
    resampler->SetReferenceImage(fixedImageReader->GetOutput());
    resampler->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType, double>::New());
    resampler->Update();

    // Write the transformed label image
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputTransformedLabelFileName);
    writer->SetInput(resampler->GetOutput());

    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject& err)
    {
        std::cerr << "ExceptionObject caught: " << err << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}