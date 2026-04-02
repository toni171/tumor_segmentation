# Liver Tumor Segmentation

A MATLAB implementation of a tumor segmenter using **Otsu's thresholding** technique, **level set evolution** and **local directional ternary pattern** technique.

## Getting Started

### Installing

In order to install the program on your local environment, open the bash inside the directory where you want to clone the repository and retrieve it using `git clone`.

```
git clone git@github.com:toni171/tumor_segmentation.git
```

### Executing program

In your MATLAB instance, navigate to the `tumor_segmentation` folder and execute the `main.m` program through the command window.

```
main
```

The program will automatically start the segmentation of all images in the dataset and save the computed metrics in the `results.json` file.

## Author

Antonio Di Leo, [@toni171](https://github.com/toni171)

## Acknowledgments

- ["Modified Otsu thresholding based level set and local directional ternary pattern technique for liver tumor segmentation"](https://link.springer.com/article/10.1007/s13198-022-01637-x), Deepak S. Uplaonkar, Virupakshappa & Nagabhushan Patil 