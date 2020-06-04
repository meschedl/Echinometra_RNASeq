#postgres table and column 
#date:05/31/2016 Parul Kudtarkar
-- Table: gene_info_table

-- DROP TABLE gene_info_table;

CREATE TABLE gene_info_table
(
  glean_id character varying(100) NOT NULL,
  spu_id character varying(100) NOT NULL,
  annotator character varying(200) NOT NULL DEFAULT 'none'::character varying,
  gene_model_check character varying(50) NOT NULL DEFAULT 'none'::character varying,
  additional_evidence character varying(50) NOT NULL DEFAULT 'none'::character varying,
  family_member text NOT NULL DEFAULT 'none'::text,
  common_name character varying(100) NOT NULL DEFAULT 'none'::character varying,
  synonyms text NOT NULL DEFAULT 'none'::text,
  best_genebank_hit character varying(100) DEFAULT 'none'::character varying,
  ortholog_homolog text NOT NULL DEFAULT 'none'::text,
  group_coordinator character varying(100) NOT NULL DEFAULT 'none'::character varying,
  difficult_annotation_categories character varying(50) NOT NULL DEFAULT 'none'::character varying,
  old_common_name character varying(200) DEFAULT 'none'::character varying,
  note character varying(150) DEFAULT 'none'::character varying,
  manual_anno character varying(5) DEFAULT 'none'::character varying,
  annotation_type character varying(50) DEFAULT 'none'::character varying,
  overlap_gene character varying(20) DEFAULT 'none'::character varying,
  whl_id character varying(100) DEFAULT 'none'::character varying,
  curator character varying(200),
  CONSTRAINT gene_info_table_spu_id_pkey PRIMARY KEY (spu_id)
)
WITH (
  OIDS=FALSE
);
ALTER TABLE gene_info_table
  OWNER TO postgres;
-- Table: comment_table

-- DROP TABLE comment_table;

CREATE TABLE comment_table
(
  glean_id character varying(100) NOT NULL,
  spu_id character varying(100) NOT NULL,
  old_comment text DEFAULT 'none'::text,
  new_comment text DEFAULT 'none'::text
)
WITH (
  OIDS=FALSE
);
ALTER TABLE comment_table
  OWNER TO postgres;
-- Table: expression_table

-- DROP TABLE expression_table;

CREATE TABLE expression_table
(
  spu_id character varying(100) NOT NULL DEFAULT 'none'::character varying,
  "time" character varying(50) NOT NULL DEFAULT 'none'::character varying,
  domain character varying(200) NOT NULL DEFAULT 'none'::character varying,
  CONSTRAINT spu_id_time_gene_domain PRIMARY KEY (spu_id, "time", domain),
  CONSTRAINT expression_table_spu_id_fkey FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE expression_table
  OWNER TO postgres;
-- Table: feature_table

-- DROP TABLE feature_table;

CREATE TABLE feature_table
(
  scaffold character varying(50) NOT NULL DEFAULT 'none'::character varying,
  source character varying(50) NOT NULL DEFAULT 'none'::character varying,
  feature character varying(50) NOT NULL DEFAULT 'none'::character varying,
  gene_start integer NOT NULL DEFAULT (-1),
  gene_end integer NOT NULL DEFAULT (-1),
  score character varying(50) NOT NULL DEFAULT 'none'::character varying,
  strand character varying(1) NOT NULL DEFAULT 'none'::character varying,
  frame character varying(2) NOT NULL DEFAULT 'none'::character varying,
  spu_id character varying(20) NOT NULL DEFAULT 'none'::character varying,
  CONSTRAINT spu_id_3_1_feature_gene_start_gene_end PRIMARY KEY (spu_id, feature, gene_start, gene_end),
  CONSTRAINT feature_table_spu_id_fkey FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE feature_table
  OWNER TO postgres;
-- Table: functional_category_table

-- DROP TABLE functional_category_table;

CREATE TABLE functional_category_table
(
  spu_id character varying(100) NOT NULL,
  class character varying(100) NOT NULL,
  sub_class character varying(100) NOT NULL,
  anno_type character varying(50) NOT NULL,
  CONSTRAINT functional_category_table_spu_id_class_sub_class_pkey PRIMARY KEY (spu_id, class, sub_class),
  CONSTRAINT functional_category_table_spu_id_fkey1 FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE functional_category_table
  OWNER TO postgres;
-- Table: go_table

-- DROP TABLE go_table;

CREATE TABLE go_table
(
  spu_id character varying(100) NOT NULL,
  ot_type character varying(100) NOT NULL,
  ot_data text NOT NULL DEFAULT 'none'::text,
  go_num character(7) NOT NULL,
  evidence_code character(3) NOT NULL DEFAULT 'none'::bpchar,
  CONSTRAINT go_table_spu_id_go_num_pkey PRIMARY KEY (spu_id, go_num),
  CONSTRAINT go_table_spu_id_fkey1 FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE go_table
  OWNER TO postgres;

-- Table: h_den_timecourse_table

-- DROP TABLE h_den_timecourse_table;

CREATE TABLE h_den_timecourse_table
(
  spu_id character varying(100) NOT NULL,
  sp_name text NOT NULL,
  hr_0 integer NOT NULL,
  hr_1 integer NOT NULL,
  hr_2 integer NOT NULL,
  hr_3 integer NOT NULL,
  hr_4 integer NOT NULL,
  hr_5 integer NOT NULL,
  hr_6 integer NOT NULL,
  hr_7 integer NOT NULL,
  hr_8 integer NOT NULL,
  hr_9 integer NOT NULL,
  hr_10 integer NOT NULL,
  hr_11 integer NOT NULL,
  hr_12 integer NOT NULL,
  hr_13 integer NOT NULL,
  hr_14 integer NOT NULL,
  hr_15 integer NOT NULL,
  hr_16 integer NOT NULL,
  hr_17 integer NOT NULL,
  hr_18 integer NOT NULL,
  hr_19 integer NOT NULL,
  hr_20 integer NOT NULL,
  hr_21 integer NOT NULL,
  hr_22 integer NOT NULL,
  hr_23 integer NOT NULL,
  hr_24 integer NOT NULL,
  hr_25 integer NOT NULL,
  hr_26 integer NOT NULL,
  hr_27 integer NOT NULL,
  hr_28 integer NOT NULL,
  hr_29 integer NOT NULL,
  hr_30 integer NOT NULL,
  hr_31 integer NOT NULL,
  hr_32 integer NOT NULL,
  hr_33 integer NOT NULL,
  hr_34 integer NOT NULL,
  hr_35 integer NOT NULL,
  hr_36 integer NOT NULL,
  hr_37 integer NOT NULL,
  hr_38 integer NOT NULL,
  hr_39 integer NOT NULL,
  hr_40 integer NOT NULL,
  hr_41 integer NOT NULL,
  hr_42 integer NOT NULL,
  hr_43 integer NOT NULL,
  hr_44 integer NOT NULL,
  hr_45 integer NOT NULL,
  hr_46 integer NOT NULL,
  hr_47 integer NOT NULL,
  hr_48 integer NOT NULL,
  CONSTRAINT h_den_timecourse_table_spu_id_sp_name_pkey PRIMARY KEY (spu_id, sp_name),
  CONSTRAINT h_den_timecourse_table_spu_id_fkey FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE h_den_timecourse_table
  OWNER TO postgres;
-- Table: input_user_table

-- DROP TABLE input_user_table;

CREATE TABLE input_user_table
(
  firstname character varying(50) NOT NULL,
  middlename character varying(50),
  lastname character varying(100) NOT NULL,
  institute character varying(150) NOT NULL,
  pi_fname character varying(50) NOT NULL,
  pi_lname character varying(100) NOT NULL,
  phone character varying(20) NOT NULL,
  email character varying(100) NOT NULL,
  pwd character varying(10) NOT NULL,
  CONSTRAINT input_user_table_email_pkey PRIMARY KEY (email)
)
WITH (
  OIDS=FALSE
);
ALTER TABLE input_user_table
  OWNER TO postgres;
-- Table: qpcr_timecourse_table

-- DROP TABLE qpcr_timecourse_table;

CREATE TABLE qpcr_timecourse_table
(
  spu_id character varying(100) NOT NULL,
  timepoint numeric(10,0) NOT NULL,
  timeunit character varying(6) NOT NULL,
  express_level numeric(10,0) DEFAULT (-1),
  CONSTRAINT qpcr_timecourse_table_spu_id_timepoint_pkey PRIMARY KEY (spu_id, timepoint),
  CONSTRAINT qpcr_timecourse_table_spu_id_fkey FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE qpcr_timecourse_table
  OWNER TO postgres;
-- Table: reagent_table

-- DROP TABLE reagent_table;

CREATE TABLE reagent_table
(
  spu_id character varying(100) NOT NULL,
  data_type character varying(30) NOT NULL,
  data_version integer NOT NULL DEFAULT 1,
  data character varying(200) NOT NULL,
  reference character varying(50) NOT NULL,
  note character varying(150) DEFAULT 'none'::character varying,
  CONSTRAINT reagent_table_spu_id_data_type_data_version_pkey PRIMARY KEY (spu_id, data_type, data_version),
  CONSTRAINT reagent_table_spu_id_fkey1 FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE reagent_table
  OWNER TO postgres;
-- Table: reference_table

-- DROP TABLE reference_table;

CREATE TABLE reference_table
(
  spu_id character varying(100) NOT NULL,
  pmid integer NOT NULL DEFAULT (-1),
  author text DEFAULT 'none'::text,
  volume character varying(30) NOT NULL DEFAULT 'none'::character varying,
  title text DEFAULT 'none'::text,
  year character varying(10) NOT NULL DEFAULT 'none'::character varying,
  journal text DEFAULT 'none'::text,
  CONSTRAINT reference_table_spu_id_pmid_pkey PRIMARY KEY (spu_id, pmid)
)
WITH (
  OIDS=FALSE
);
ALTER TABLE reference_table
  OWNER TO postgres;
-- Table: seq_table

-- DROP TABLE seq_table;

CREATE TABLE seq_table
(
  spu_id character varying(20) NOT NULL,
  version character varying(3) NOT NULL,
  method character varying(50) NOT NULL,
  info_type character varying(10) NOT NULL,
  info text NOT NULL DEFAULT 'none'::text,
  sub_version character varying(8) NOT NULL,
  CONSTRAINT seq_table_pkey PRIMARY KEY (spu_id, version, info_type, sub_version),
  CONSTRAINT seq_table_spu_id_fkey FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE seq_table
  OWNER TO postgres;
-- Table: spatial_expression_table

-- DROP TABLE spatial_expression_table;

CREATE TABLE spatial_expression_table
(
  spu_id character varying(100) NOT NULL,
  space_time character varying(20) NOT NULL,
  expression character varying(20) NOT NULL,
  confidence character varying(20) NOT NULL,
  CONSTRAINT spacial_expression_table_spu_id_space_time_pkey PRIMARY KEY (spu_id, space_time),
  CONSTRAINT spacial_expression_table_spu_id_fkey FOREIGN KEY (spu_id)
      REFERENCES gene_info_table (spu_id) MATCH SIMPLE
      ON UPDATE NO ACTION ON DELETE NO ACTION
)
WITH (
  OIDS=FALSE
);
ALTER TABLE spatial_expression_table
  OWNER TO postgres;
-- Table: whl_table

-- DROP TABLE whl_table;

CREATE TABLE whl_table
(
  spu_id character varying(100) NOT NULL,
  whl_id character varying(100) NOT NULL DEFAULT 'none'::character varying,
  CONSTRAINT gene_info_table_spu_whl_id_pkey PRIMARY KEY (spu_id, whl_id)
)
WITH (
  OIDS=FALSE
);
ALTER TABLE whl_table
  OWNER TO postgres;
