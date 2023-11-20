**YYC for file version control system** is a DNA storage codec algorithm suitable for encoding files in different 
versions. It is adjusted based on YYC. It can transcode two adjacent binary sequences into one DNA sequence, and only 
recode the modified parts of each version file. This algorithm can help to achieve a high-density, high-feasibility DNA 
storage of multi-version files based on DNA synthesis.

## Environment Configuration
The kit is developed by **Python3.9**.

In addition, the packages we are calling now is as follows:

- [x] sys
- [x] os
- [x] random
- [x] math
- [x] copy
- [x] struct
- [x] datetime
- [x] pickle

## Kit Tree Diagram
```html
├── examples                          // Test module
│    ├── files                        // Test files
│    │    ├── match_files             // Match files
│    │    │    ├── version2           // Match files for the second version file
│    │    │    ├── version3           // Match files for the third version file
│    │    │    ├── version4           // Match files for the fourth version file
│    │    ├── version1.txt            // Original version file in txt format
│    ├── output                       // Generated files from encoding or decoding
├── yyc
│    ├── utils                        // Util module
│    │    ├── data_handle.py          // Conversion of DNA motifs and binary document
│    │    ├── index_operator.py       // Processing the relationship between index and data
│    │    ├── log.py                  // Output the logs in console
│    │    ├── model_saver.py          // Save model to file and load model from file
│    │    ├── monitor.py              // Get the progress situation and the time left
│    │    ├── validity.py             // Determining whether a DNA sequence is easy or not for sequencing and synthesis
│    ├── pipeline.py                  // Main calling function
│    ├── scheme.py                    // YYC for file version control system
├── README.md                         // Description document of kit
```