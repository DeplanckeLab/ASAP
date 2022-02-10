package model;

public enum Mode
{
	Preparsing,
	Parsing,
	PreparseMetadata,
	ParseMetadata,
	RegenerateNewOrganism,
	CreateCellSelection,
	CreateGODB,
	CreateKeggDB,
	DifferentialExpression,
	FindMarkers,
	MarkerEnrichment,
	Normalization,
	Scaling,
	UpdateEnsemblDB,
	Enrichment,
	ModuleScore,
	DimensionReduction,
	IndexByCell,
	GetIndex,
	GetGeneStats,
	ExtractRow,
	ExtractCol,
	ExtractDataset,
	ListMetadata,
	ExtractMetadata,
	MatchValues,
	RemoveMetaData,
	CopyMetaData,
	FilterCols,
	FilterRows,
	FilterDEMetadata;
	
	public static String toArrayString()
	{
		StringBuffer sb = new StringBuffer("[");
		String prefix = "";
		for(Mode m:Mode.values())
		{
			sb.append(prefix).append(m.toString());
			prefix = ", ";
		}
		sb.append("]");
		return sb.toString();
	}
	


}
