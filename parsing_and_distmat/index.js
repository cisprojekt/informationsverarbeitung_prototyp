//Global variable to store the data from the csv file
let hammingData = Array();
//Parsing part of the code
//
//
// Define a class to store the data points
class DataPoint {
    data = "";
    flag = "";
    constructor(data, flag) {
        this.data = data;
        this.flag = flag;
    }    
}

function handleFile() {
    // Get the input element and the selected file
    var fileInput = document.getElementById('fileInput');
    var file = fileInput.files[0];

    if (file) {
        var reader = new FileReader();

        // Set up the onload callback function
        reader.onload = function (e) {
            // Get the contents of the file
            var contents = e.target.result;

            // Split the CSV content into lines considering CR and LF as line endings
            var lines = contents.split(/\r?\n/);

            // Assuming the first line contains headers
            var headers = lines[0].split(',');

            // Initialize an array to store the objects
            var dataPointArray = [];

            // Iterate over the remaining lines
            for (var i = 1; i < lines.length; i++) {
                var values = lines[i].split(',');
                console.log(values);
                var dataPoint = new DataPoint(values[0], values[1]);
                dataPointArray.push(dataPoint);
            }
            console.log(dataPointArray); 
            console.log(dataPointArray[0].data); // Print the field 'data'
            console.log("done");
            hammingData = dataPointArray.map(point => point.data);
            console.log(hammingData);
        };
        // Read the file as text
        reader.readAsText(file);
    } else {
        alert('Please select a file.');
    }
}

//ensure WASM is loaded
//this is an alternative to the Module.onRuntimeInitialized = function(){}; method
let wasmReady = new Promise(resolve => {
    Module.onRuntimeInitialized = resolve;
});


async function onCalcButtonClick() {
    await wasmReady;
    console.log("Euclidean test")
    const matrix = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    const array = matrix.flat();
    console.log(matrix.length);
    console.log(matrix[0].length);
    //console.log(array);
    
    let buf = Module._malloc(array.length * Float32Array.BYTES_PER_ELEMENT);

    Module.HEAPF32.set(array, buf / Float32Array.BYTES_PER_ELEMENT);

    //call the distance matrix function
    //returns the pointer to the result array
    let resultPtr = Module.ccall('calculateEuclideanDistanceMatrix', 'number', ['number','number', 'number'], [buf, matrix.length, matrix[0].length]);

    //create a typed array from the pointer containing the distmat as flattened array
    let distmat = new Float32Array(Module.HEAPF32.buffer, resultPtr, array.length);

    for (let i = 0; i < (matrix.length*(matrix.length+1)/2); i++) {
        console.log(distmat[i]);
    }

    Module._free(buf);

    console.log("Hamming test")
    const nucleotides = hammingData;
    const string_array = nucleotides.map(str => new Uint8Array(str.split('').map(c => c.charCodeAt(0))));

    for (let i = 0; i < string_array.length; i++) {
        console.log(string_array[i]);
    }

    // Allocate memory for the strings and copy the character codes to the allocated memory
    const charPtrs = string_array.map(chars => {
        const ptr = Module._malloc(chars.length * chars.BYTES_PER_ELEMENT);
        Module.HEAPU8.set(chars, ptr);
        return ptr;
    });

    // Allocate memory for the array of pointers
    const ptrBuf = Module._malloc(charPtrs.length * Int32Array.BYTES_PER_ELEMENT);

    // Copy the array of pointers to the allocated memory
    Module.HEAP32.set(charPtrs, ptrBuf / Int32Array.BYTES_PER_ELEMENT);

    //call the distance matrix function
    //returns the pointer to the result array
    let resultPtr2 = Module.ccall('calculateHammingDistanceMatrix', 'number', ['number','number','number'], [ptrBuf, nucleotides.length, nucleotides[0].length]);

    //create a typed array from the pointer containing the distmat as flattened array
    let hamdistmat = new Int32Array(Module.HEAP32.buffer, resultPtr2, nucleotides.length*(nucleotides.length+1)/2);
    
    for (let i = 0; i < (nucleotides.length*(nucleotides.length+1)/2); i++) {
        console.log(hamdistmat[i]);
    }

    Module._free(ptrBuf);
    
    //Display the matrices
    //
    //
    //

    function displayMatrix(matrix, container) {
    // Create a table
    const table = document.createElement('table');

    // For each row in the matrix
    for (let i = 0; i < matrix.length; i++) {
        // Create a table row
        const row = document.createElement('tr');

        // For each element in the row
        for (let j = 0; j < matrix[i].length; j++) {
            // Create a table cell
            const cell = document.createElement('td');

            // Set the cell text to the matrix element
            cell.textContent = matrix[i][j];

            // Add the cell to the row
            row.appendChild(cell);
        }

        // Add the row to the table
        table.appendChild(row);
    }

    // Add the table to the container
    container.appendChild(table);
    }   
    matrix2 = new Array(nucleotides.length);
        for (let i = 0; i < nucleotides.length; i++) {
            matrix2[i] = new Array(nucleotides.length);
        }
        let  k = 0;
        for (let i = 0; i < nucleotides.length; i++) {
            for (let j = 0; j < nucleotides.length; j++) {
                if(j < i + 1) {
                    matrix2[i][j] = hamdistmat[k];
                    k++;
                }
                else matrix2[i][j] = 0;
            }
        }
    
    const container = document.getElementById('matrix');
    displayMatrix(matrix2, container);
    
};