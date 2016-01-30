# Ludwig
Looking for negative selection in cancer genomes - DPhil rotation project with Benjamin Schuster-Böckler at the Ludwig Institute for Cancer Reasearch at the University of Oxford

In theory the driver script dNdS.sh should automatically download the raw data, analyse it and output the results:
```
$ dNdS.sh
```

## Data Flow Diagram
![Data Flow Diagram](data_flow_diagram.png)

## Dependencies
Python 2.6.6 (r266:84292)
R 3.2.2
	data.table 1.9.6
	ggplot2 2.0.0
	metRology 0.9-17
	ggrepel 0.4
graphviz 2.38.0 (20140413.2041) (dot for creating data flow diagram)