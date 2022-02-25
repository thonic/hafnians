pipeline {
    agent any
    stages {
        stage('Run hafnian.py') {
            steps {
                sh '''
                    pwd
                    conda env list
                '''
            }
        }
    }
}