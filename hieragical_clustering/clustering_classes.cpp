#include <iostream>
#include <vector>

class DistanceMatrix {

	private:
		unsigned int dim;
		std::vector<std:vector<double> > datapoints; //Lieber Eigen::Matrix ?	
		//obere Dreiecksmatrix um Redundanz zu vermeiden
		std:vector<double> datapoints_alt; //eindimensionale Repräsentation mit Indexmanagment

	public:
		DistamceMatrix();
		~DistanceMatrix();

}

DistanceMatrix::DistanceMatrix() { //Distanzmatrix aus Stream also hochgeladenem File
	this->data = convert_stream_to_distmat();
	dim = data.size();

}

DistanceMatrix::DistanceMatrix(double **xy) { //Distanzmatrix aus xy-Wertetabelle

	this->data = convert_xy_table_to_distmat();
	dim = data.size();

DistanceMatrix::~DistanceMatrix() {
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


}

Cluster::Cluster(unsigned int ID) {
	this->ID = ID;
}

Cluster::~Cluster() {

}

class Dendrogram {

}

Dendrogram calc_dendrogram(zoomlevel, rawdata, distance_update, algorithm) {
	//Implementieren der Cython-Funktion(en)
}

ClusterMap extract_set_from_cut_dendrogram(zoomlevel, dendrogram);
	//irgendeine Variante, wie man das dendrogramm entsprechend der Zoomstufe an der passenden Stelle abschneidet
	//und die übrig bleibenden Mengen als Vektor von Clustertypen zurückgibt
