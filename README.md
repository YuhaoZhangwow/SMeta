
# SMeta Tool

This command-line tool facilitates the processing of FASTA sequences based on the provided metadata and sequence files. It's designed to be flexible, allowing various settings through optional parameters to meet different processing needs.

## Getting Started

These instructions will guide you through the setup process and demonstrate how to run `fasta.exe` on your local machine for development and testing purposes.

### Prerequisites

- Ensure you have a compatible Windows environment to run `.exe` files.
- Prepare the metadata and sequence files in formats that `fasta.exe` can process.

### Installation

Download the latest release of `fasta.exe` from the Releases section of this GitHub repository and save it to a suitable directory on your computer.

## Usage

To run the tool, you'll need to open a command prompt and navigate to the directory where `fasta.exe` is located. The output will be generated in result.txt

### Basic Command

```
fasta.exe -meta <metadata_path> -single <sequence_path> [optional parameters]
```

### Parameters

#### Mandatory Parameters
- **-meta/-m`<string>`**: Path to the metadata file [Mandatory].
- **-single/-s`<string>`**: Path to the single sequence directory, the program will read all files under this path [Mandatory].

#### Optional Parameters
- **-help/-h**: Help information.
- **-slice `<number>`**: Number of segment trees per single cell sequence you would like to generate for processing (default: `2`).
- **-maxedges `<number>`**: Maximum number of edges built for one metagenome sequence during alignment inside metagenome sequences(default: `200`).
- **-maxedgeseg `<number>`**: Maximum number of edges per metagenome sequence built during alignment with single cell sequences and metagenome sequences(default: `50`).
- **-debug/-d`<1 or 0>`**: Debug mode, enable with `1` and disable with `0`, `1` will provide progress bar and other information like the number of edges built (default: `1`).
- **-maxp `<double>`**: Maximum threshold for the proportion of points that Label Propagation Algorithm should cover(default: `0.95`).
- **-time/-t/-T `<1 or 0>`**: Enable timing in process (default: `0`).

### Examples

**Running with Required Parameters Only:**

```bash
fasta.exe -meta "C:/data/metadata.txt" -single "C:/data/sequence.fasta"
```

**Running with All Parameters:**

```bash
fasta.exe -meta "C:/data/metadata.txt" -single "C:/data/sequence.fasta" -slice 3 -maxedges 250 -maxedgeseg 75 -debug 1 -maxp 0.92 -t 1
```

### Compilation For environment other than Windows

To compile the `fasta.cpp` source file into an executable, follow the instructions below. These steps assume you are using `g++` on a Windows system with paths to necessary libraries.

#### Prerequisites

Ensure you have `g++` installed and accessible in your environment. Additionally, the following libraries should be installed and properly configured:
- Eigen (3.4.0)
- Boost (1.84.0)

You may need to adjust the include paths depending on where these libraries are installed on your system.

#### Compile Command

Open your command prompt and run the following command:

```bash
g++.exe <path to cpp file> -o <path where you want to output> -g3 -std=c++17 -O3 -fopenmp -I"PATH_TO\MinGW64\include" -I"PATH_TO\eigen-3.4.0\eigen-3.4.0" -I"PATH_TO\boost_1_84_0\boost_1_84_0" -L"PATH_TO\MinGW64\lib" -static-libgcc -g3
```

#### Explanation of Compile Flags

- `-o`: Specifies the output file name for the compiled executable.
- `-g3`: Enables debug information to be generated at the highest level to facilitate debugging.
- `-std=c++17`: Specifies the C++ standard to C++17 for using modern language features.
- `-O3`: Enables compiler optimizations for faster execution.
- `-fopenmp`: Enables the OpenMP for parallel programming.
- `-I`: Specifies the include directories for header files.
- `-L`: Specifies the directories where libraries are searched.
- `-static-libgcc`: Links the static version of the GCC lib, allowing the executable to run on systems without the GCC installed.

## Contributing

Contributions to this project are welcome! Please consider the following ways to contribute:

1. Reporting bugs or issues.
2. Suggesting enhancements or new features.
3. Contributing to the code via pull requests.

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests to us.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

- Grateful that I am still alive
