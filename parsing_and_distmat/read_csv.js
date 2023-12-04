class DataPoint {
    data = "";
    flag = "";
    constructor(data, flag) {
        this.data = data;
        this.flag = flag;
    }
    calculateDistance(dataPoint) {
        return Math.sqrt(this.data**2 - dataPoint.data**2);
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
        };
        // Read the file as text
        reader.readAsText(file);
    } else {
        alert('Please select a file.');
    }
}