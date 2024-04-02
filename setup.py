from setuptools import setup, find_packages


setup(
    name="localhgt",
    version="1.0.0",
    author="wangshuai",
    author_email="wshuai294@gmail.com",
    description="Detect HGT from microbiome.",
    url="https://github.com/deepomicslab/LocalHGT", 
    packages=find_packages(),


    scripts=['scripts/build_UHGG_reference.py', 'scripts/infer_HGT_breakpoint.py','scripts/infer_HGT_event.py']
)
