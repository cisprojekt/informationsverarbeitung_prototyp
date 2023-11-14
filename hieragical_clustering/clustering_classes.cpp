#include <iostream>
#include <vector>


//
class DistanceMatrix {

	private:
		unsigned int dim;
		double *datapoints; //Lieber Eigen::Matrix ?	
		//obere Dreiecksmatrix um Redundanz zu vermeiden
		// std:vector<double> datapoints_alt; eindimensionale Repräsentation mit Indexmanagment

	public:
		DistamceMatrix();
		~DistanceMatrix();
		double *getmatrix();

}

DistanceMatrix::DistanceMatrix() { //Distanzmatrix aus Stream also hochgeladenem File
	this->data = convert_stream_to_distmat();
	dim = data.size();

}

DistanceMatrix::DistanceMatrix(double *xy, unsigned int n) { //Distanzmatrix aus xy-Wertetabelle, xy hat folgende Struktur: Array bestehend pointern mit 2 double Werten [[2.43, 4.32], [1.23, 7.11], ...

	this->datapoints = convert_xy_table_to_distmat(xy, n);
	this->dim = n;

DistanceMatrix::~DistanceMatrix() {
}

double *convert_xy_table_to_distmat(double *xy_table, unsigned int n) {

	double *distmat = new double[(n*(n-1)/2];
	int k, i, j;
	for (i=k=0; i<n, i++) {
		for (j=i+1; j<n, j++) {
			distmat[k] = sqrt(pow((xy_table[i][0]+xy_table[j][0]), 2) + pow((xy_table[i][1]+xy_table[j][1])                     //bis hierhin ist die Funktion auch ganz gut als generische Funktion -> könnte ausgelagert werden. Noch generischer mit distance(x[i], x[j]) stattt sqrt ...
			k++;
		}
	}
	return distmat;
}

double *getmatrix() {
	return this->datapoints;
}			

class ClusterMap {

	private:
		unsigned char zoomlevel; //Idee ist die Erzeugung einer ClusterMap Instanz pro Zoomlevel
		// damit kann man in Zukunft dynamisch entweder nur eine Instanz berechnen, oder alle Zoomlevel vorberechnen
		unsigned char mode; //auf welche Art wurde das clustering durchgeführt, mode 1 = hierarchisch_agglomerativ
		DistanceMatrix rawdata; //enthält die Distanzmatrix der Ursprungsdatenpunkte
		DistanceMatrix cluster_dist; //enthält die Distanzen zwischen den erstellten Clustern;
	public:
		
		ClusterMap(ClusterMap cm1);
		ClusterMap(unsigned char mode);
		std::vector<Cluster> global; //beinhaltet alle (berechneten) Cluster aus den Eingabedaten
		void converge_cluster_to_xy(); //für die map-Schnittstelle
		DistanceMatrix calc_distmat(); //berechnet Distanzmatrix der fertigen Map

		
}

ClusterMap::ClusterMap(ClusterMap cm1) {
	//Berechnung der Clustermap aus einer anderen (bspw. mit höherer Zoomstufe)
	//Modus der Berechnung ergibt sich aus cm1->mode
}
ClusterMap::ClusterMap(unsigned char mode, unsigned char submode, unsigned char subsubmode) {
	//Berechnung der Clustermap aus der zugrundeliegenden Distanzmatrix
	if (mode == 1) {
		Dendrogram dendrogram = calc_dendrogram(zoomlevel, rawdata, distance_update, algorithm) //im wesentlichen die Cython Funktion 
	    global = extract_sets_from_cut_dendrogram(zoomlevel, dendrogram); //dendrogramm an der richtigen Stelle abschneiden und die Mengendaten die darin gespeichert werden in einen Vektor aus Cluster-Datentyp verwandeln
		cluster_dist = extract_distances_from_cut_dendrogram(dendrogram); 
		//aus dem abgeschnittenen Dendrogramm die Distanzen in eine neue Distanzmatrix schreiben

	}
}
class Cluster {
	//noch nicht ganz sicher wofür wir die Klasse genau brauchen. Momentan denke ich, dass sie vor allem sinnvoll
	// ist, um später Metadaten zu speichern, die wir ja im fortgeschrittenen Stadium der Entwicklung anzeigen lassen wollen
	private:
		unsigned int ID;
		unsigned int num_elements;
		//std::vector<Cluster> elements; 
			//eventuell naive Variante alle Clusterinhalte zu speichern
			// Problem: mit dem zehnten Clusteringlevel ist das ein zehnfach gewrappter Vektor - super ineffizient
			

	public:
		Cluster(unsigned int ID);
		~Cluster();
		unsigned int getID();

}

Cluster::Cluster(unsigned int ID) {
	this->ID = ID;
}

Cluster::~Cluster() {

}

unsigned int getID() {
			return this->ID;
		}

class Dendrogram {  //bzw. linkage-matrix

}

Dendrogram calc_dendrogram(zoomlevel, rawdata, distance_update, algorithm) {
	//Implementieren der Cython-Funktion(en)
}

ClusterMap extract_set_from_cut_dendrogram(zoomlevel, dendrogram);
	//irgendeine Variante, wie man das dendrogramm entsprechend der Zoomstufe an der passenden Stelle abschneidet
	//und die übrig bleibenden Mengen als Vektor von Clustertypen zurückgibt
