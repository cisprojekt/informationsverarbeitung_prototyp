/* global Module */
/* global msgpack5 */

// maximum tries for module loading
const MAX_TRIES = 10;
let numTries = 0;

// set waiting interval
const moduleInterval = setInterval(() => {
	
	// increment number of tries of module not loaded
    if (!Module) {
        numTries++;
    }

    if (numTries >= MAX_TRIES) {
        clearInterval(moduleInterval);
    }

	// module is loaded
    if (Module && Module.calledRun && msgpack5) {
        clearInterval(moduleInterval);

		// load messagepack
        const msgpack = msgpack5();

		// creates random points
		// with x and y coordinate (int)
        function create_rand_points(n) {
            const points = [];
            for (let i = 0; i < n; i++) {
                const point = {
                    positionId: "position" + (i + 1),
                    x: parseInt(Math.random() * 500),
                    y: parseInt(Math.random() * 500)
                };
                points.push(point)
            }
            return points;
        }

		// create 40 random points
		nPoints = 40;
        var myPoints = create_rand_points(nPoints);
		
		console.log(myPoints);

		// allocate memory in the heap for the encoded data
        const encodedPositions = msgpack.encode(myPoints);
        let inputPointer = Module._malloc(encodedPositions.length);
        Module.HEAP8.set(
            encodedPositions,
            inputPointer / encodedPositions.BYTES_PER_ELEMENT
        );

		// call the C++ function to cluster the points
        const outputPointer = Module._malloc(8);
        const processedPointsIdsPointer = Module.ccall(
            "processPoint",
            "number",
            ["number", "number", "number", "number"],
            [inputPointer, encodedPositions.length, outputPointer, nPoints]
        );

        let offset = Module.getValue(outputPointer, "i64");
		
		// convert BigInt to number
		offset = Number(offset)
		
		// create subarray
		// TODO: endpointer not working for n > 40
        const processedPointsIdsData = new Uint8Array(
            Module.HEAPU8.subarray(
                processedPointsIdsPointer,
                processedPointsIdsPointer + offset
            )
        );

		// decode using messagepack
        const processedPointsIds = msgpack.decode(processedPointsIdsData);
		
		console.log(processedPointsIds);
		
		// print the points and their cluster id
        var content = "";
		for (let i = 0; i < processedPointsIds.length; i++) {
			content += "Point " + i + " is in cluster " + processedPointsIds[i]+ "<br>";
		}
        document.getElementById("result").innerHTML = content;

		// deallocate memory
        Module._free(inputPointer);
        Module._free(outputPointer);
        Module._free(processedPointsIdsPointer);
    }
}, 10);