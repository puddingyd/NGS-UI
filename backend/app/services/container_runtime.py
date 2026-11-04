"""Keep container bindings local to the current managed job."""
from pathlib import Path

BIND_VARS = ('SINGULARITY_BINDPATH', 'SINGULARITY_BIND', 'APPTAINER_BINDPATH', 'APPTAINER_BIND')


def job_bind_environment(environment, work_root, required):
    """Preserve reference binds, replace temp destinations, discard job staging.

    Tertiary children can outlive their parent's disposable runtime directory.
    Never carry that directory into an independent PRS/Nextflow invocation,
    even while it still exists. Deduplicate by container destination rather
    than the complete source:destination string.
    """
    root = Path(work_root).absolute()
    transient = [root/'contexts'/f'shared-{mode}'/'.runtime/jobs' for mode in ('dragen', 'nckuh')]
    bindings = {}
    for variable in BIND_VARS:
        for value in environment.get(variable, '').split(','):
            value = value.strip()
            if not value:
                continue
            pieces = value.split(':')
            source = Path(pieces[0]).absolute()
            destination = pieces[1] if len(pieces) > 1 and pieces[1] else pieces[0]
            if destination in ('/tmp', '/var/tmp') or any(source.is_relative_to(p) for p in transient):
                continue
            bindings[destination] = value
    for value in required:
        pieces = value.split(':')
        destination = pieces[1] if len(pieces) > 1 and pieces[1] else pieces[0]
        bindings[destination] = value
    combined = ','.join(bindings.values())
    # *_BIND and *_BINDPATH can both be consumed by Apptainer. Use a single
    # binding variable per runtime; keep the legacy alias for Singularity.
    return dict(APPTAINER_BIND='', SINGULARITY_BIND='',
                APPTAINER_BINDPATH=combined, SINGULARITY_BINDPATH=combined)


def without_bind_destinations(environment, destinations):
    """The probe passes these exact mounts on argv; avoid adding them twice."""
    result = dict(environment)
    for name in BIND_VARS:
        values = []
        for value in result.get(name, '').split(','):
            if not value: continue
            pieces = value.split(':')
            target = pieces[1] if len(pieces) > 1 and pieces[1] else pieces[0]
            if target not in destinations: values.append(value)
        result[name] = ','.join(values)
    return result
