#include "maestro.h"

#include <klib/rc.h>

#include <vdb/manager.h>
#include <vdb/schema.h>
#include <sra/sraschema.h> /* VDBManagerMakeSRASchema */
#include <vdb/database.h>
#include <vdb/table.h>
#include <kdb/table.h>
#include <vdb/vdb-priv.h> /* VDBManagerGetKDBManagerRead */
#include <kdb/meta.h> /* KMetadata */

using namespace std;

// Return zero if we can't access the SRA metadata for any reason. 
uint64_t number_of_bases(const string &m_accession)
{
	uint64_t ret = 0;

	// This code is based on the SRA toolkit program sra-stat.c
	const VDBManager* vmgr = NULL;
	const VDatabase* db = NULL;
	const VTable* vtbl = NULL;
	VSchema* schema = NULL;
	const KTable* ktbl = NULL;
	const KMetadata* meta = NULL;
	const KMDataNode* node = NULL;
	const KMDataNode* subnode = NULL;

	try{
		if(VDBManagerMakeRead(&vmgr, NULL) != 0){
			throw __FILE__ ":number_of_bases: VDBManagerMakeRead != 0";
		}

		if(VDBManagerMakeSRASchema(vmgr, &schema) != 0){
			throw __FILE__ ":number_of_bases: VDBManagerMakeSRASchema != 0";
		}

		rc_t rc = VDBManagerOpenTableRead( vmgr, &vtbl, schema, "%s", m_accession.c_str() );

		if( (rc != 0) && (GetRCObject(rc) == (enum RCObject)rcTable) &&
        	(GetRCState(rc) == rcIncorrect)){
			
			// SRA is not a table
			if(VDBManagerOpenDBRead( vmgr, &db, schema, m_accession.c_str() ) != 0){
				throw __FILE__ ":number_of_bases: VDBManagerOpenDBRead != 0";
			}

			if(VDatabaseOpenTableRead(db, &vtbl, "%s", "SEQUENCE") != 0){
				throw __FILE__ ":number_of_bases: VDatabaseOpenTableRead != 0";
			}

			rc = 0;
		}

		if(rc != 0){
			throw __FILE__ ":number_of_bases: Unable to open table";
		}

		if(VTableOpenKTableRead(vtbl, &ktbl) != 0){
			throw __FILE__ ":number_of_bases: VTableOpenKTableRead != 0";
		}

		if(KTableOpenMetadataRead(ktbl, &meta) != 0){
			throw __FILE__ ":number_of_bases: KTableOpenMetadataRead != 0";
		}

		if(KMetadataOpenNodeRead(meta, &node, "%s", "STATS/TABLE") != 0){
			throw __FILE__ ":number_of_bases: KMetadataOpenNodeRead != 0";
		}

		if(KMDataNodeOpenNodeRead(node, &subnode, "%s", "BASE_COUNT") != 0){
			throw __FILE__ ":number_of_bases: KMDataNodeOpenNodeRead != 0";
		}

		if(KMDataNodeReadAsU64(subnode, &ret) != 0){
			throw __FILE__ ":number_of_bases: KMDataNodeReadAsU64 != 0";
		}

		// Not currently using the BIO_BASE_COUNT table
		//readStatsMetaNode(node, "STATS/TABLE", "BIO_BASE_COUNT", &ret, quick, false);
	}
	catch(...){
		ret = 0;
	}

	if(subnode != NULL){
		KMDataNodeRelease(subnode);
	}

	if(node != NULL){
		KMDataNodeRelease(node);
	}

	if(meta != NULL){
		KMetadataRelease(meta);
	}

	if(ktbl != NULL){
		KTableRelease(ktbl);
	}

	// Note that we need to release vtbl, not tbl
	if(vtbl != NULL){
		VTableRelease(vtbl);
	}

	if(db != NULL){
		VDatabaseRelease(db);
	}

	if(schema != NULL){
		VSchemaRelease(schema);
	}

	if(vmgr != NULL){
		VDBManagerRelease(vmgr);
	}

	return ret;
}