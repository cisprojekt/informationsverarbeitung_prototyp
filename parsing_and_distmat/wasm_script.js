
let wasmExports = null;

let wasmMemory = new WebAssembly.Memory({initial: 65.536, maximum: 65.536});

let wasmTable = new WebAssembly.Table({
    'initial': 1,
    'maximum': 1,
    'element': 'anyfunc'
});

let asmLibraryArg = { 
    "__handle_stack_overflow": ()=>{},
    "emscripten_resize_heap": ()=>{},
    "__lock": ()=>{}, 
    "__unlock": ()=>{},
    "memory": wasmMemory, 
    "table": wasmTable 
};

var info = {
    'env': asmLibraryArg,
    'wasi_snapshot_preview1': asmLibraryArg
};

async function loadWasm(){
    let response = await fetch('distmat.wasm');
    let bytes = await response.arrayBuffer();
    let wasmObj = await WebAssembly.instantiate(bytes, info);
    wasmExports = wasmObj.instance.exports;
}

loadWasm();

/*
// Fetch and instantiate the WebAssembly module
fetch('distmat.wasm')
.then(response => response.arrayBuffer())
.then(bytes => WebAssembly.instantiate(bytes))
.then(results => {
    // Get the WebAssembly instance
    const instance = results.instance;

    // Your matrix
    const matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];

    // Flatten the matrix into a 1D array
    const array = matrix.flat();

    // Allocate memory in the WebAssembly module's memory space
    const ptr = instance.exports._malloc(array.length * 4); // Assume each element is 4 bytes (e.g., int32)

    // Create a new Uint8Array backed by the WebAssembly memory
    const wasmMemory = new Uint8Array(instance.exports.memory.buffer);

    // Copy the array data into the allocated memory
    for (let i = 0; i < array.length; i++) {
        const baseAddress = ptr + i * 4; // Assume each element is 4 bytes
        new DataView(wasmMemory.buffer, baseAddress, 4).setFloat32(0, array[i], true); // Assume float32 data
    }

    // Call the WebAssembly function
    instance.exports.calculateDistanceMatrix(ptr, matrix.length, matrix[0].length);

    // Free the allocated memory (if the WebAssembly module has a _free function)
    if (instance.exports._free) {
        instance.exports._free(ptr);
    }
});
*/
/*
function testfunc(memory){
    console.log("hellooooooooooooooooooooo")
    const matrix = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    const array = matrix.flat();
    const array2 = new Float64Array(memory.buffer, 0, array.length*4)
    array2.set(array)
    var ptr = wasmExports.calculateDistanceMatrix(array2.byteOffset, 3, 3);
    console.log(ptr);}*/
