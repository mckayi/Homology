(* ::Package:: *)

SetDirectory["C:\\Users\\Ian\\Desktop\\Persistent Homology DNA\\Data"]

fileNames = FileNames["*.bcalm.fa"]

filehandle = StringTake[fileNames,10]

For[i=1, i<=Length[fileNames],i++,
(
	SetDirectory[StringJoin["C:\\Users\\Ian\\Desktop\\Persistent Homology DNA\\Output\\",filehandle[[i]]]];
	h=Import[StringJoin[Directory[],"\\",fileNames[[i]],".mma.txt"]];
	g=ToExpression[h];
	gsub = Subgraph[g, Select[ConnectedComponents[g], Length[#] == Max[Length /@ ConnectedComponents[g]] &]];
	hd3=GraphEmbedding[gsub,"HighDimensionalEmbedding",3];
	hd10=GraphEmbedding[gsub,"HighDimensionalEmbedding",10];
	hd31=GraphEmbedding[gsub,"HighDimensionalEmbedding",31];
	Export[StringJoin[filehandle[[i]],".hd3.mat"],hd3];
	Export[StringJoin[filehandle[[i]],".hd10.mat"],hd10];
	Export[StringJoin[filehandle[[i]],".hd31.mat"],hd31];
)
]

Exit[]
