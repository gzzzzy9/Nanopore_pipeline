from setuptools import setup, find_packages

setup(
    name=nanopore4repseq,
    version=1.0.0,
    packages=find_packages(where=src),
    package_dir={ src},
    include_package_data=True,
    install_requires=[
        pandas,
        numpy,
        streamlit,
        plotly,
    ],
    entry_points={
        'console_scripts' [
            # 格式：命令名 = 包名.模块名函数名
            'nanopore-run = nanopore4repseq.mainmain',
            'nanopore-ui = nanopore4repseq.appmain', 
        ],
    },
)