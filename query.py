from typing import Optional, Union
from astroquery.gaia import Gaia
table = Gaia.load_table('gaiadr3.gaia_source')

def freeQuery(data: Optional[Union[str,list]], conditions: Optional[Union[str,list]]):
    dat=''
    if isinstance(data, list):
        dat = ''
        for i in range(len(data)-1):
            dat += data[i]+', '
        dat += data[len(data)-1]
    else: dat=data

    if isinstance(conditions, str):
        conditions = conditions.split(',')

    base_query = """SELECT 
    {data}
    FROM gaiadr3.gaia_source
    WHERE {cond}
    """

    cond = ''
    for i in range(len(conditions)-1):
        cond += conditions[i]+"""
        AND """
    cond += conditions[len(conditions)-1]
    query = base_query.format(data=dat, cond=cond)

    job = Gaia.launch_job(query)
    res = job.get_results()

    return res
#########################################################################################
si può fare un self._basequery = base_query nell '__init__'
e da li fare tutte le funzioni di query
