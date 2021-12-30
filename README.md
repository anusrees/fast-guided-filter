# Fast Guided filter

  
This code is based on the guided filter implementation using OpenCV as given here: https://github.com/atilimcetin/guided-filter. This is a pure C++ based implementation and uses multithreading for faster execution speed.
Guided filter is an edge-preserving smoothing filter like the bilateral filter. It is straightforward to implement and has linear complexity independent of the kernel size. For more details about this filter see [[Kaiming10]](http://research.microsoft.com/en-us/um/people/kahe/eccv10/).

  
  

## Usage

  

The interface consists of one simple function `guidedFilterMono` and a class `AuxMono` for single channel, and function and class `guidedFilterColor` and `AuxColor` respectively. If you have multiple images to filter with the same guidance image then use `AuxMono` or `AuxColor` class to avoid extra computations on initialization stage. The code supports single-channel and 3-channel (color) guidance images and `uint8_t`, `int8_t`, `uint16_t`, `int16_t`, `uint32_t`,`int32_t`, `float` and `double` data types.

  
  

## Example

  

These examples are adapted from the [original MATLAB implementation](http://research.microsoft.com/en-us/um/people/kahe/eccv10/guided-filter-code-v1.rar).

  

  
  

### Flash/no-flash denoising

  

```c++
int  r = 8, numThreads=atoi(argv[1]);
double  eps = 0.02 * 0.02 * 255 * 255;

cv::Mat  I = cv::imread("../cave-noflash.png", cv::IMREAD_COLOR);
cv::Mat  p = cv::imread("../cave-flash.png", cv::IMREAD_COLOR);

cv::Mat  q(p.cols, p.rows, p.type());

AuxColor<uint8_t> *ac = new  AuxColor<uint8_t>(I.cols, I.rows, p.channels(), I.channels(), numThreads);

auto  t1 = std::chrono::high_resolution_clock::now();

guidedFilterColor<uint8_t>(I.data, p.data, r, eps, q.data, ac);

auto  t2 = std::chrono::high_resolution_clock::now();

std::cout  <<  std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() <<  " ms\n";

cv::imwrite("../output.png", q);

delete  ac;
```

  

[![Cave Flash](http://atilimcetin.com/guided-filter/img_flash/cave-flash-small.png)](http://atilimcetin.com/guided-filter/img_flash/cave-flash.png)

[![Cave No Flash](http://atilimcetin.com/guided-filter/img_flash/cave-noflash-small.png)](http://atilimcetin.com/guided-filter/img_flash/cave-noflash.png)

[![Cave Denoised](http://atilimcetin.com/guided-filter/img_flash/cave-denoised-small.png)](http://atilimcetin.com/guided-filter/img_flash/cave-denoised.png)

  
  
Check main.cpp for more details.

### Building the program
Instructions for running the program on Linux
1) Install cmake, g++ and OpenCV (Has dependency on OpenCV for reading and writing the image files only):
`sudo apt-get install cmake g++ libopencv-dev`
2) Using CmakeList.txt to build executable:
`mkdir build && cd build && cmake -S .. && make`
3) For running the program, execute the following command from the build folder:
`.\main <num_of_threads>`

## License

  

MIT license.