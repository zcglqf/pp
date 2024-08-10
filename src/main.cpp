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
            << " <AtlasImage> <TargetImage> <AtlasLabel> <OutputTransformedLabel>"
            << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Press any key to continue:";
    std::cin.get();

    const char* atlasImageFileName = argv[1];
    const char* targetImageFileName = argv[2];
    const char* atlasLabelFileName = argv[3];
    const char* outputTransformedLabelFileName = argv[4];

    // Read the atlas (fixed) and target (moving) images
    ReaderType::Pointer atlasImageReader = ReaderType::New();
    atlasImageReader->SetFileName(atlasImageFileName);
    atlasImageReader->Update();

    ReaderType::Pointer targetImageReader = ReaderType::New();
    targetImageReader->SetFileName(targetImageFileName);
    targetImageReader->Update();

    ReaderType::Pointer atlasLabelReader = ReaderType::New();
    atlasLabelReader->SetFileName(atlasLabelFileName);
    atlasLabelReader->Update();

    // Set up the registration
    RegistrationType::Pointer registration = RegistrationType::New();
    MetricType::Pointer metric = MetricType::New();
    OptimizerType::Pointer optimizer = OptimizerType::New();
    TransformType::Pointer transform = TransformType::New();

    registration->SetMetric(metric);
    registration->SetOptimizer(optimizer);
    registration->SetFixedImage(atlasImageReader->GetOutput());  // Atlas image as fixed
    registration->SetMovingImage(targetImageReader->GetOutput());  // Target image as moving
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

    // Resample the atlas label image to match the target image space
    using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(atlasLabelReader->GetOutput());
    resampler->SetTransform(finalTransform);
    resampler->SetUseReferenceImage(true);
    resampler->SetReferenceImage(atlasImageReader->GetOutput());
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