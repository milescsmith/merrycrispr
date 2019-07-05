import json
import sys
import time
import os
import gzip
import shutil
from distutils.spawn import find_executable
from multiprocessing import cpu_count
from subprocess import check_call
from tempfile import mkdtemp, mkstemp
from typing import Optional, Dict, Tuple
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import urlopen, Request

import requests
import tqdm

"""Class for retreiving information from ensembl.org
Heavily based on the example given at
`https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client`
"""


class EnsemblRestClient:
    def __init__(self, server: str = "http://rest.ensembl.org", reqs_per_sec: int = 15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(
        self,
        endpoint: str,
        hdrs: Optional[Dict[str, str]] = None,
        params: Optional[Dict[str, str]] = None,
    ):
        if hdrs is None:
            hdrs = {}

        if "Content-Type" not in hdrs:
            hdrs["Content-Type"] = "application/json"

        if params:
            endpoint += "?" + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if "Retry-After" in e.headers:
                    retry = e.headers["Retry-After"]
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    "Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n".format(
                        endpoint, e
                    )
                )
        return data

    """Retrieve information for a species from Ensembl
    \f
    Parameters
    -----------
    self : :class:`species_getter.EnsemblRestClient`
        Calling object
    value : `str`
        Attribute selection value.  For instance, if selecting by "name" (i.e. scientific name) 
        then "bos_taurus" or if selecting by "accession" then "GCA_002263795.2"
    attribute : `str`, optional (default: 'name')
        Attribute to use when selecting the species. Acceptable values include 'taxon_id', 
        'accession', 'aliases', 'division', 'groups', 'release', 'name', 'strain', 
        'strain_collection', 'display_name', 'assembly', and 'common_name'

    Return
    -------
    `dict`
    """

    def get_species_info(
        self, value: str, attribute: Optional[str] = "name"
    ) -> Dict[str, str]:
        all_species = self.perform_rest_action(endpoint="/info/species")
        potential_matches = [_ for _ in all_species["species"] if value in _[attribute]]
        if not potential_matches:
            raise ValueError(
                f"No species with a {attribute} matching {value} was found"
            )
        elif len(potential_matches) > 1:
            print(
                f"Multiple values where {attribute} matches {value} were found, "
                f"so you will need to be more specific.  The following matches were found:"
            )
            for index, item in enumerate(potential_matches):
                print(index)
                for _ in item:
                    print(f"{_}: {item[_]}")
            sys.exit(1)
        else:
            return potential_matches[0]

    """Generalized getter to retrieve files from ftp.ensembl.org
    \f
    Parameters
    -----------
    self : :class:`species_getter.EnsemblRestClient`
        Calling object
    file_url : `str`
        Location of the remote file.
    destination : `str`, default=False
        Folder in which to place the downloaded file.
    force_redownload : `bool`, default=False

    Return
    ------
    `str` 
        Folder in which the file was placed
    """

    def get_genome_file(
        self,
        file_url: str,
        destination: Optional[str] = None,
        force_redownload: bool = False,
    ) -> str:
        if not destination:
            destination = mkstemp()

        try:
            r = requests.get(file_url, stream=True)
        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)

        file_size = int(r.headers["Content-Length"])
        if (
            os.path.exists(destination)
            and os.path.isfile(destination)
            and os.stat(destination).st_size == file_size
            and not force_redownload
        ):
            print(f"File is already present and at {destination}")
        else:
            chunk_size = 1024
            num_bars = int(file_size / chunk_size)
            with open(destination, "wb") as fp:
                for chunk in tqdm.tqdm(
                    r.iter_content(chunk_size=chunk_size),
                    total=num_bars,
                    unit="KB",
                    desc=destination,
                    leave=True,
                ):
                    fp.write(chunk)
        return destination

    """Retrieve the latest GTF for a given species
    \f
    Parameters
    -----------
    self : :class:`species_getter.EnsemblRestClient`
        Calling object
    species_value : `str`
        Value being used to select species of interest, i.e. "canis_familiaris", "CanFam3.1", 
        "9615", "GCA_000002285.2", etc...
    species_attribute : `str`
        Value being used to select species of interest, i.e. name, assembly, taxon_id, accession, 
        etc...
    destination : `str`, optional.
        Folder in which to place the downloaded file.  If none provided,
        the file will be downloaded to a temp directory.

    Return
    -------
    `str` 
        Folder in which the file was placed
    """

    def get_annotation(
        self,
        species_value: str,
        species_attribute: str,
        destination: Optional[str] = None,
    ) -> str:
        species_info = self.get_species_info(
            value=species_value, attribute=species_attribute
        )
        destination = destination.rstrip("/")
        annotation_url = (
            f"https://ftp.ensembl.org/pub/release-{species_info['release']}"
            f"/gtf/{species_info['name']}"
            f"/{species_info['name'].capitalize()}."
            f"{species_info['assembly']}."
            f"{species_info['release']}.gtf.gz"
        )
        return self.get_genome_file(
            file_url=annotation_url,
            destination=f"{destination}/"
            f"{species_info['name'].capitalize()}."
            f"{species_info['assembly']}."
            f"{species_info['release']}.gtf.gz",
        )

    """Retrieve the latest primary DNA sequence FASTA for a given species
    \f
    Parameters
    -----------
    self : :class:`species_getter.EnsemblRestClient`
        Calling object
    species : `str`
        Species of interest.  Either scientific or common name.
    destination : `str`, optional.
        Folder in which to place the downloaded file.  If none provided,
        the file will be downloaded to a temp directory.
    
    Return
    -------
    `str` 
        Folder in which the file was placed
    """

    def get_sequences(
        self,
        species_value: str,
        species_attribute: str,
        destination: Optional[str] = None,
    ) -> str:
        species_info = self.get_species_info(
            value=species_value, attribute=species_attribute
        )
        destination = destination.rstrip("/")
        annotation_url = (
            f"https://ftp.ensembl.org/pub/release-{species_info['release']}"
            f"/fasta/{species_info['name']}"
            f"/dna/{species_info['name'].capitalize()}."
            f"{species_info['assembly']}.dna.primary_assembly.fa.gz"
        )
        return self.get_genome_file(
            file_url=annotation_url,
            destination=f"{destination}/"
            f"{species_info['name'].capitalize()}."
            f"{species_info['assembly']}.dna.primary_assembly.fa.gz",
        )


"""List the species currently available from Ensembl
\f
Parameters
-----------

Return
-------
`None`
"""


def available_species() -> None:
    client = EnsemblRestClient()
    species = client.perform_rest_action(endpoint="/info/species")
    if species:
        for _ in species["species"]:
            print(f"{_['common_name']}: {_['name']}")


"""Retrieve the latest GTF for a given species
\f
Parameters
-----------
species : `str`
    Species of interest.  Either common or scientific name.
resource_folder : `str`, optional (default: None)
    Folder to download files to.  If not provided, a new directory called "resource_folder"
    will be created and used.

Return
-------
:class:`typing.Tuple`[`str`,`str`]
    The locations of the files downloaded.
"""


def get_resources(
    species_value: str, species_attribute: str, resource_folder: Optional[str] = None
) -> Tuple[str, str]:
    client = EnsemblRestClient()
    resource_folder = os.path.expanduser(resource_folder)
    if not os.path.exists(resource_folder):
        try:
            os.makedirs("resource_folder")
        except:
            ValueError(f"problem making {resource_folder}")

    gtf = client.get_annotation(
        species_value=species_value,
        species_attribute=species_attribute,
        destination=resource_folder,
    )
    fasta = client.get_sequences(
        species_value=species_value,
        species_attribute=species_attribute,
        destination=resource_folder,
    )
    return gtf, fasta


"""Build an index for Bowtie
\f
Parameters
-----------
fasta : `str`
    Sequence file (in FASTA format) to build index against
dest: `str`, optional (default: None).
    Directory in which to place index.
cpus: `int`, optional (default: 0).
    Number of CPUs to use in building index.

Return
-------
`str` 
    Folder in which the file was placed
"""


def build_bowtie_index(fasta: str, dest: Optional[str] = None, cpus: int = 0) -> str:
    if cpus == 0:
        cpus = cpu_count()
    program = find_executable("bowtie-build")

    # Bowtie cannot index a gzipped fasta, which is how they come from Ensembl
    if fasta[-3:] == ".gz":
        with gzip.open(fasta, 'rb') as fgz:
            with open(fasta[:-3], "wb") as fa:
                shutil.copyfileobj(fgz, fa)
        fasta = fasta[:-3] # set the name to the uncompressed item

    if not dest:
        dest = mkdtemp()
    command = f"{program} -f --threads {cpus} {fasta} {dest}"

    try:
        check_call(command.split())
    except BaseException:
        raise SystemExit("Bowtie encountered an error. Please check the log file.")
    return dest
