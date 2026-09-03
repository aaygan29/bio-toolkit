from __future__ import annotations

import os
import sqlite3
from dataclasses import dataclass, field


@dataclass
class ProteinEntry:
    accession: str
    gene_name: str
    full_name: str
    sequence: str
    organism: str = "Homo sapiens"
    isoforms: list[dict] = field(default_factory=list)
    variants: list[dict] = field(default_factory=list)
    domains: list[dict] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.sequence)


@dataclass
class VariantEntry:
    accession: str
    variant_id: str
    position: int
    ref_aa: str
    alt_aa: str
    description: str = ""
    pathogenicity: str = "unknown"


DIAGNOSTIC_KMER_SIZE = 9

BUILTIN_PROTEINS: dict[str, dict] = {
    "P04637": {
        "gene": "TP53",
        "name": "Cellular tumor antigen p53",
        "sequence": (
            "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP"
            "DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK"
            "SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE"
            "RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS"
            "SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
            "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG"
            "GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
        ),
        "isoforms": [
            {"id": "P04637-2", "name": "Isoform 2 (p53 beta)"},
            {"id": "P04637-3", "name": "Isoform 3 (p53 gamma)"},
            {"id": "P04637-4", "name": "Isoform 4 (delta N p53)"},
        ],
        "variants": [
            {"id": "rs1042522", "pos": 72, "ref": "R", "alt": "P",
             "desc": "Pro72Arg polymorphism", "path": "benign"},
            {"id": "rs28934571", "pos": 175, "ref": "R", "alt": "H",
             "desc": "R175H Li-Fraumeni", "path": "pathogenic"},
            {"id": "rs28934574", "pos": 248, "ref": "R", "alt": "W",
             "desc": "R248W common mutation", "path": "pathogenic"},
        ],
        "domains": [
            {"name": "Transactivation", "start": 1, "end": 61},
            {"name": "Proline-rich", "start": 64, "end": 92},
            {"name": "DNA-binding", "start": 102, "end": 292},
            {"name": "Tetramerization", "start": 324, "end": 356},
        ],
    },
    "P00533": {
        "gene": "EGFR",
        "name": "Epidermal growth factor receptor",
        "sequence": (
            "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEV"
            "VLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALA"
            "VLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDF"
            "QNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGC"
            "TGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYV"
            "VTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFK"
            "NCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAF"
            "ENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKL"
            "FGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCN"
            "LLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVM"
            "GENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVV"
            "ALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGS"
            "GAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGI"
            "CLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAA"
            "RNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSY"
            "GVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPK"
            "FRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQ"
            "GFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDS"
            "IDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNT"
            "VQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVA"
            "PSSSEFIGA"
        ),
        "isoforms": [
            {"id": "P00533-2", "name": "Isoform 2 (EGFRvIII)"},
            {"id": "P00533-3", "name": "Isoform 3 (secreted)"},
        ],
        "variants": [
            {"id": "rs121913444", "pos": 858, "ref": "L", "alt": "R",
             "desc": "L858R activating (NSCLC)", "path": "pathogenic"},
            {"id": "rs121913465", "pos": 719, "ref": "G", "alt": "S",
             "desc": "G719S activating", "path": "pathogenic"},
            {"id": "rs121913437", "pos": 790, "ref": "T", "alt": "M",
             "desc": "T790M resistance", "path": "pathogenic"},
        ],
        "domains": [
            {"name": "Receptor L domain", "start": 57, "end": 168},
            {"name": "Furin-like", "start": 176, "end": 341},
            {"name": "Tyrosine kinase", "start": 712, "end": 979},
        ],
    },
    "P38398": {
        "gene": "BRCA1",
        "name": "Breast cancer type 1 susceptibility protein",
        "sequence": (
            "MDLSALRVEEVQNVINAMQKILECPICLELLIKEPVSTKCDHIFCKFCMLKLLNQKKGPS"
            "QCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLEYANSYNFAKKENNSPEHLK"
            "DEVSKIQSMVRNVMQHCKEGIFKKVPQSNHVQESSSEESMSVREQPLNTEDEGSDSVPSP"
            "SLQSALRSEAPVRTSSENPCHERQSSREPPGSAETGQASSENASDASASPEPKELSSSTE"
            "PKHKKPLKTPAQKRESKKEPLQINGPAIPESPLAPKENAKNQTERLQESLFLKRKEEEEQ"
            "TVQKNNQTPPESASLPENLEDAEVKQNFLKLPENKHLREDKEKELQSSVMEKNSKLFESK"
            "PKSSSGSGESSSEGSLSSESESLVLLDASEPQNLPSSCSNDLDPLVYSSVIQDSINSQEL"
            "SNTDEIDGLAANSNLQTHSPSSINKHSPTLTIDLNSSVKSLEVNGDLSREILNSSPSMLC"
            "TPTEKCSSKEGEEPDIQNQKNKKFPNISRDSQITNQGLNALSTCSNSDHFQFTATNKSDP"
            "SASFQWRNSGSSEQDLDFSGNSNISLENLEQLNPLGIDRGSLSKILSEPENLKVNGLSEN"
            "QNTRPESLVVSDLGSPEDNLVSKLQNPSNNSLLILGSEEGVTVEKPVMKMFSKPMHSINP"
            "NKGLKPAFSMSVLSECSDPSKILEEVNEPSEAEYVDPSLEENDLTSSQKVDASDKHPQNE"
            "TSSPSVKAKALNQSGNSSAEEVLLSSSSPTSQVNQSQLPLKSSGVSSQKVPCGVKEYHIE"
            "DPESKANQLVMDTEEILVKSEETLVTSSNKNSKSTEQETQTTTNKISVLPSNSKEEHEKS"
            "LTQKSNTQREELPSKNSKEDITASNVDILSSKECMKADGQVDNSSKKAESIRNTLTASKA"
            "KNKENCSTSKTEEKNFGNNNDENLSVNSSNKNDVKSGKNSKRESNSKGPEKESAVNTSEE"
            "KDSKTHKEFSHISSLGSKAEQEDRNSNNAETGSENENKEASESPSSNTIKAAEEPRNNSS"
            "SKEFIPSSREHKKELQKADLSNDSHKKFSVANINELGEKSNTEINAEERKGKISKANPDK"
            "SYKSVVVDLNFDKKSKLNEGHQAYVDTATTEESCFTTEKSRRPDSHSSHQSQKEIKSANS"
            "ASSRTLSNNNSEKDPGHKDSFKDTPASSSNIEQLTDQENSSQMSDSVQKDSLVNGKNNSS"
            "VSVCTSGSSIGSSGFSSSSSSFSSESGSGSSGSESSSESSSESSSSLSSESSGSSSESEE"
            "PSVSSPLTRKGSSVEASTKLMTEPASCRIRLKEQDSSKNSIRRIHTTGGKMNSSERHVDV"
            "RKKNMTLSRSGINQTAETGGHSKAEGKSCSSILSNSSIYDRNSSSGDKEMFESVRSSKSL"
            "FNSDLTSVPLNIGINKKEASSSWQKKESFMNNNEMLLSRGQQTSKNKEEQSSPSKNDQNT"
            "LTDQQQLEKRESISKIDDYDQKAFTEKARESSSKKCLMVMSLQSEKQLSSCSSETLNRQV"
            "SSVSSGSSVSNIPAPISMPWEDSEGSSVTNSQEVNVFKKVRGRRLRK"
        ),
        "isoforms": [
            {"id": "P38398-2", "name": "Isoform 2 (BRCA1a)"},
            {"id": "P38398-3", "name": "Isoform 3 (BRCA1b)"},
        ],
        "variants": [
            {"id": "rs80357906", "pos": 1775, "ref": "Q", "alt": "*",
             "desc": "Q1775* (185delAG founder)", "path": "pathogenic"},
        ],
        "domains": [
            {"name": "RING finger", "start": 24, "end": 64},
            {"name": "BRCT 1", "start": 1642, "end": 1731},
            {"name": "BRCT 2", "start": 1756, "end": 1847},
        ],
    },
    "P68871": {
        "gene": "HBB",
        "name": "Hemoglobin subunit beta",
        "sequence": (
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSNADAVMGNPK"
            "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG"
            "KEFTPPVQAAYQKVVAGVANALAHKYH"
        ),
        "isoforms": [],
        "variants": [
            {"id": "rs334", "pos": 6, "ref": "E", "alt": "V",
             "desc": "E6V sickle cell", "path": "pathogenic"},
            {"id": "rs33969966", "pos": 26, "ref": "G", "alt": "D",
             "desc": "HbE G26D", "path": "pathogenic"},
        ],
        "domains": [{"name": "Globin", "start": 1, "end": 147}],
    },
    "P69905": {
        "gene": "HBA1",
        "name": "Hemoglobin subunit alpha",
        "sequence": (
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK"
            "KVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAV"
            "HASLDKFLAVSTVLTSKYR"
        ),
        "isoforms": [],
        "variants": [
            {"id": "rs414120", "pos": 94, "ref": "D", "alt": "N",
             "desc": "HbA2' (D94N)", "path": "benign"},
        ],
        "domains": [{"name": "Globin", "start": 1, "end": 142}],
    },
    "P02768": {
        "gene": "ALB",
        "name": "Serum albumin",
        "sequence": (
            "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFE"
            "DHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPER"
            "NECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAK"
            "RYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLS"
            "QRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEK"
            "PLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYS"
            "VVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKF"
            "QNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLH"
            "EKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKK"
            "QTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"
        ),
        "isoforms": [],
        "variants": [],
        "domains": [{"name": "Albumin domain 1", "start": 19, "end": 114},
                     {"name": "Albumin domain 2", "start": 191, "end": 292},
                     {"name": "Albumin domain 3", "start": 382, "end": 486}],
    },
    "P60709": {
        "gene": "ACTB",
        "name": "Actin, cytoplasmic 1",
        "sequence": (
            "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQS"
            "KRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPTLLTEAPLNPKANREKMT"
            "QIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDL"
            "AGRDLTDYLMKILTERGYSFVTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSY"
            "ELPDGQVITIGNERFRCPETLFQPSFIGMESAGIHETTYNSIMKCDIDIRKDLYANNVLS"
            "GGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISKQ"
            "EYDESGPSIVHRKCF"
        ),
        "isoforms": [],
        "variants": [],
        "domains": [{"name": "Actin", "start": 1, "end": 375}],
    },
    "P01308": {
        "gene": "INS",
        "name": "Insulin",
        "sequence": (
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED"
            "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
        ),
        "isoforms": [],
        "variants": [],
        "domains": [
            {"name": "Insulin B chain", "start": 90, "end": 110},
            {"name": "Insulin A chain", "start": 25, "end": 54},
        ],
    },
    "P42212": {
        "gene": "GFP",
        "name": "Green fluorescent protein",
        "organism": "Aequorea victoria",
        "sequence": (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
            "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV"
            "NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLAD"
            "HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        ),
        "isoforms": [],
        "variants": [
            {"id": "S65T", "pos": 65, "ref": "S", "alt": "T", "desc": "Enhanced brightness", "path": "benign"},
        ],
        "domains": [{"name": "GFP beta-barrel", "start": 1, "end": 238}],
    },
    "P13569": {
        "gene": "CFTR",
        "name": "CF transmembrane conductance regulator",
        "sequence": (
            "MQRSPLEKASVVSKLFFSWTRPILRKGYRQRLELSDIYQIPSVDSADNLSEKLEREWDRE"
            "LASKKNPKLINALRRCFFWRFMFYGIFLYLGEVTKAVQPLLLGRIIASYDPDNKEERSIA"
            "IYLGIGLCLLFIVRTLLLHPAIFGLHHIGMQMRIAMFSLIYKKTLKLSSRVLDKISIGQL"
            "VSLLSNNLNKFDEGLALAHFVWIAPLQVALLMGLIWELLQASAFCGLGFLIVLALFQAGL"
            "GRMMMKYRDQRAGKISERLVITSEMIENIQSVKAYCWEEAMEKMIENLRQTELKLTRKAA"
            "YVRYFNSSAFFFSGFFVVFLSVLPYALIKGIILRKIFTTISFCIVLRMAVTRQFPWAVQH"
            "WYDSLGAINKIQDFLQKQEYKTLEYNLTTTEVVMENVTAFWEEGFGELFEKAKQNNNNRK"
            "TSNGDDSLFFSNFSLLGTPVLKDINFKIERGQLLAVAGSTGAGKTSLLMMIMGELEPSEG"
            "KIKHSGRISFCSQFSWIMPGTIKENIIFGVSYDEYRYRSVIKACQLEEDISKFAEKDNIV"
            "LGEGGITLSGGQRARISLARAVYKDADLYLLDSPFGYLDVLTEKEIFESCVCKLMANKTR"
            "ILVTSKMEHLKKADKILILHEGSSYFYGTFSELQNLQPDFSSKLMGCDSFDQFSAERRNS"
            "ILTETLHRFSLEGDAPVSWTETKKQSFKQTGEFGEKRKNSILNPINSIRKFSIVQKTPLQ"
            "MNGIEEDSDEPLERRLSLVPDSEQGEAILPRISVISTGPTLQARRRQSVLNLMTHSVNQG"
            "QNIHRKTTASTRKVSLAPQANLTELDIYSRRLSQETGLEISEEINEEDLKECFFDDMESI"
            "PAVTTWNTYLRYITVHKSLIFVLIWCLVIFLAEVAASLVVLWLLGNTPLQDKGNSTHSRN"
            "NSYAVIITSTSSYYVFYIYVGVADTLLAMGFFRGLPLVHTLITVSKILHHKMLHSVLQAP"
            "MSTLNTLKAGGILNRFSKDIAILDDLLPLTIFDFIQLLLIVIGAIAVVAVLQPYIFVATV"
            "PVIVAFIMLRAYFLQTSQQLKQLESEGRSPIFTHLVTSLKGLWTLRAFGRQPYFETLFHK"
            "ALNLHTANWFLYLSTLRWFQMRIEMIFVIFFIAVTFISILTTGEGEGRVGILTLAMNIMS"
            "TLQWAVNSSIDVDSLMRSVSRVFKFIDMPTEGKPTKSTKPYKNGQLSKVMIIENSHVKKD"
            "DIWPSGGQMTVKDLTAKYTEGGNAILENISFSISPGQRVGLLGRTGSGKSTLLSAFLRLL"
            "NTEGEIQIDGVSWDSITLQQWRKAFGVIPQKVFIFSGTFRKNLDPYEQWSDQEIWKVADE"
            "VGLRSVIEQFPGKLDFVLVDGGCVLSHGHKQLMCLARSVLSKAKILLLDEPSAHLDPVTY"
            "QIIRRTLKQAFADCTVILCEHRIEAMLECQQFLVIEENKVRQYDSIQKLLNERSLFRQAI"
            "SPSDRVKLFPHRNSSKCKSKPQIAALKETTEEEVQDTRL"
        ),
        "isoforms": [],
        "variants": [
            {"id": "rs113993960", "pos": 508, "ref": "F", "alt": "del",
             "desc": "F508del major CF", "path": "pathogenic"},
        ],
        "domains": [
            {"name": "NBD1", "start": 433, "end": 585},
            {"name": "R domain", "start": 595, "end": 860},
            {"name": "NBD2", "start": 1225, "end": 1415},
        ],
    },
    "P04406": {
        "gene": "GAPDH",
        "name": "Glyceraldehyde-3-phosphate dehydrogenase",
        "sequence": (
            "MGKVKVGVNGFGRIGRLVTRAAFNSGKVDIVAINDPFIDLNYMVYMFQYDSTHGKFHGTV"
            "KAENGKLVINGNPITIFQERDPSKIKWGDAGAEYVVESTGVFTTMEKAGAHLQGGAKRVI"
            "ISAPSADAPMFVMGVNHEKYDNSLKIISNASCTTNCLAPLAKVIHDNFGIVEGLMTTVHA"
            "ITATQKTVDGPSGKLWRDGRGALQNIIPASTGAAKAVGKVIPELNGKLTGMAFRVPTANV"
            "SVVDLTCRLEKPAKYDDIKKVVKQASEGPLKGILGYTEHQVVSSDFNSDTHSSTFDAGAG"
            "IALNDHFVKLISWYDNEFGYSNRVVDLMAHMASKE"
        ),
        "isoforms": [],
        "variants": [],
        "domains": [{"name": "GAPDH", "start": 1, "end": 335}],
    },
    "P01375": {
        "gene": "TNF",
        "name": "Tumor necrosis factor",
        "sequence": (
            "MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQR"
            "EEEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVEL"
            "RDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQR"
            "ETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
        ),
        "isoforms": [],
        "variants": [],
        "domains": [{"name": "TNF domain", "start": 77, "end": 233}],
    },
}


class ProteinRepository:
    def __init__(self, db_path: str | None = None):
        if db_path is None:
            db_path = os.path.join(os.path.dirname(__file__), "protein_repo.db")
        self._conn: sqlite3.Connection | None = None
        self._db_path = db_path
        self._init_once()

    @property
    def conn(self) -> sqlite3.Connection:
        if self._conn is None:
            self._conn = sqlite3.connect(self._db_path)
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def _init_once(self) -> None:
        c = self.conn
        existing = c.execute(
            "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='proteins'"
        ).fetchone()[0]
        if existing:
            return
        c.executescript("""
            CREATE TABLE proteins (
                accession TEXT PRIMARY KEY, gene_name TEXT NOT NULL,
                full_name TEXT, organism TEXT DEFAULT 'Homo sapiens',
                sequence TEXT NOT NULL, sequence_length INTEGER
            );
            CREATE TABLE isoforms (
                id TEXT PRIMARY KEY, parent_accession TEXT NOT NULL,
                name TEXT, differences TEXT,
                FOREIGN KEY (parent_accession) REFERENCES proteins(accession)
            );
            CREATE TABLE variants (
                id TEXT, accession TEXT NOT NULL, position INTEGER NOT NULL,
                ref_aa TEXT, alt_aa TEXT, description TEXT,
                pathogenicity TEXT DEFAULT 'unknown',
                PRIMARY KEY (accession, id)
            );
            CREATE TABLE domains (
                accession TEXT NOT NULL, name TEXT NOT NULL,
                start_pos INTEGER NOT NULL, end_pos INTEGER NOT NULL,
                FOREIGN KEY (accession) REFERENCES proteins(accession)
            );
            CREATE TABLE diagnostic_kmers (
                kmer TEXT NOT NULL, accession TEXT NOT NULL,
                position INTEGER NOT NULL, is_unique INTEGER DEFAULT 0,
                PRIMARY KEY (kmer, accession)
            );
            CREATE INDEX idx_variants_pos ON variants(position);
            CREATE INDEX idx_kmers ON diagnostic_kmers(kmer);
        """)
        self._seed_builtin()
        c.commit()

    def _seed_builtin(self) -> None:
        for acc, data in BUILTIN_PROTEINS.items():
            seq = data["sequence"].replace(" ", "").replace("\n", "")
            self.conn.execute(
                "INSERT OR IGNORE INTO proteins VALUES (?,?,?,?,?,?)",
                (acc, data["gene"], data["name"],
                 data.get("organism", "Homo sapiens"), seq, len(seq)),
            )
            for iso in data.get("isoforms", []):
                self.conn.execute(
                    "INSERT OR IGNORE INTO isoforms VALUES (?,?,?,?)",
                    (iso["id"], acc, iso["name"], iso.get("differences", "")),
                )
            for var in data.get("variants", []):
                self.conn.execute(
                    "INSERT OR IGNORE INTO variants VALUES (?,?,?,?,?,?,?)",
                    (var["id"], acc, var["pos"], var["ref"], var["alt"],
                     var.get("desc", ""), var.get("path", "unknown")),
                )
            for dom in data.get("domains", []):
                self.conn.execute(
                    "INSERT OR IGNORE INTO domains VALUES (?,?,?,?)",
                    (acc, dom["name"], dom["start"], dom["end"]),
                )
            self._compute_diagnostic_kmers(acc, seq)
        self.conn.commit()

    def _compute_diagnostic_kmers(self, accession: str, seq: str) -> None:
        k = DIAGNOSTIC_KMER_SIZE
        other_seqs = {
            a: d["sequence"].replace(" ", "").replace("\n", "")
            for a, d in BUILTIN_PROTEINS.items() if a != accession
        }
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            is_unique = 0 if any(kmer in other for other in other_seqs.values()) else 1
            self.conn.execute(
                "INSERT OR IGNORE INTO diagnostic_kmers VALUES (?,?,?,?)",
                (kmer, accession, i, is_unique),
            )

    def get_by_accession(self, accession: str) -> ProteinEntry | None:
        row = self.conn.execute(
            "SELECT * FROM proteins WHERE accession = ?", (accession,)
        ).fetchone()
        if not row:
            return None
        return self._build_entry(row)

    def get_by_gene(self, gene_name: str) -> list[ProteinEntry]:
        rows = self.conn.execute(
            "SELECT * FROM proteins WHERE gene_name = ?", (gene_name.upper(),)
        ).fetchall()
        return [self._build_entry(r) for r in rows]

    def search_by_sequence(self, query: str, min_identity: float = 0.7) -> list[tuple[ProteinEntry, float, int]]:
        k = DIAGNOSTIC_KMER_SIZE
        query_kmers = {
            query[i:i+k].upper()
            for i in range(max(0, len(query) - k + 1))
        }
        if not query_kmers:
            return []

        rows = self.conn.execute(
            f"SELECT DISTINCT accession FROM diagnostic_kmers WHERE kmer IN ({','.join('?'*len(query_kmers))})",
            list(query_kmers),
        ).fetchall()

        results: list[tuple[ProteinEntry, float, int]] = []
        for row in rows:
            entry = self.get_by_accession(row["accession"])
            if not entry:
                continue
            covered = min(len(query), len(entry.sequence))
            matches = sum(1 for i in range(covered) if query[i] == entry.sequence[i])
            identity = matches / max(covered, 1)
            if identity >= min_identity:
                results.append((entry, identity, matches))

        results.sort(key=lambda x: x[1], reverse=True)
        return results

    def get_domains(self, accession: str) -> list[dict]:
        return [
            {"name": r["name"], "start": r["start_pos"], "end": r["end_pos"]}
            for r in self.conn.execute(
                "SELECT name, start_pos, end_pos FROM domains WHERE accession=? ORDER BY start_pos",
                (accession,),
            ).fetchall()
        ]

    def get_variants(self, accession: str) -> list[VariantEntry]:
        return [
            VariantEntry(
                accession=r["accession"], variant_id=r["id"],
                position=r["position"], ref_aa=r["ref_aa"], alt_aa=r["alt_aa"],
                description=r["description"], pathogenicity=r["pathogenicity"],
            )
            for r in self.conn.execute(
                "SELECT * FROM variants WHERE accession=?", (accession,)
            ).fetchall()
        ]

    def compare_isoforms(self, accession: str) -> list[dict]:
        return [
            dict(r) for r in self.conn.execute(
                "SELECT * FROM isoforms WHERE parent_accession=?", (accession,)
            ).fetchall()
        ]

    def _build_entry(self, row: sqlite3.Row) -> ProteinEntry:
        acc = row["accession"]
        return ProteinEntry(
            accession=acc,
            gene_name=row["gene_name"],
            full_name=row["full_name"],
            sequence=row["sequence"],
            organism=row["organism"],
            variants=[
                {"id": v.variant_id, "position": v.position, "ref": v.ref_aa,
                 "alt": v.alt_aa, "description": v.description, "pathogenicity": v.pathogenicity}
                for v in self.get_variants(acc)
            ],
            isoforms=[dict(r) for r in self.conn.execute(
                "SELECT * FROM isoforms WHERE parent_accession=?", (acc,)).fetchall()
            ],
            domains=self.get_domains(acc),
        )

    def list_all(self) -> list[ProteinEntry]:
        return [self._build_entry(r) for r in self.conn.execute("SELECT * FROM proteins").fetchall()]

    def close(self) -> None:
        if self._conn:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
